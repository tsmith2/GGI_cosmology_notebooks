#!/usr/bin/env python3
"""Colab-friendly wrapper for the simplified_CMB_lite TT solver.

Model parameter vector:

    theta = [A_s, omega_cdm, omega_b, H0, n_s, Delta_N_eff]

where H0 is in km/s/Mpc, and omega_cdm and omega_b are the usual physical
density parameters. Delta_N_eff is the deviation from the default N_eff=3.
The likelihood functions also accept ell_min and ell_max, with defaults
ell_min=2 and ell_max=2500.
"""

from __future__ import annotations

import atexit
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np


PACKAGE_DIR = Path(__file__).resolve().parents[1]
CPP_DIR = PACKAGE_DIR / "cpp"
EXE = CPP_DIR / "two_fluid_tt"
CACHE_DIR = PACKAGE_DIR / ".cache"
FIGURE_DIR = PACKAGE_DIR / "figures"
MOCK_DATA_FILE = PACKAGE_DIR / "data" / "mock_tt_planck_like_lite_neutrino.npz"
BINNING_FILE = PACKAGE_DIR / "data" / "planck_tt_binning_lmax2500.npz"
BESSEL_X_CACHE = CACHE_DIR / "lite_bessel_x_l2500_dx0p1_x15000.bin"

MPLCONFIG_DIR = PACKAGE_DIR / ".mplconfig"
MPLCONFIG_DIR.mkdir(exist_ok=True)
CACHE_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_DIR))

import matplotlib.pyplot as plt


PARAM_NAMES = [
    r"$A_s$",
    r"$\omega_{\rm cdm}$",
    r"$\omega_b$",
    r"$H_0$",
    r"$n_s$",
    r"$\Delta N_{\rm eff}$",
]

FID_THETA = np.array([2.1e-9, 0.1201, 0.0223, 67.0, 0.965, 0.0], dtype=float)
LCDM_BOUNDS = np.array(
    [
        [1.0e-9, 3.5e-9],
        [0.09, 0.16],
        [0.017, 0.028],
        [55.0, 80.0],
        [0.85, 1.08],
        [-2.9, 5.0],
    ],
    dtype=float,
)

ELL_MIN = 2
ELL_MAX = 2500
N_SOURCE = 1001
K_SAMPLES_PER_PERIOD = 40.0
N_THREADS = 0

_BUILT_CPP = False
_SERVERS: dict[tuple[int, int, float, int], "TwoFluidCPPServer"] = {}


@dataclass(frozen=True)
class TTSpectrum:
    ell: np.ndarray
    cl: np.ndarray
    dell_over_2pi: np.ndarray


def _validate_theta(theta: np.ndarray) -> np.ndarray:
    """Return theta as a float array after checking the H0 convention."""
    theta = np.asarray(theta, dtype=float)
    if theta.shape != (len(FID_THETA),):
        raise ValueError("theta must have six entries")
    if not np.isfinite(theta[3]) or theta[3] < 10.0:
        raise ValueError(
            "theta[3] must be H0 in km/s/Mpc, for example 67.0. "
            "Do not pass 0.67 here."
        )
    if 3.0 + theta[5] < 0.0:
        raise ValueError("N_eff = 3 + Delta_N_eff must be non-negative")
    return theta.copy()


def ensure_cpp_executable(*, force_rebuild: bool = False) -> None:
    """Compile the C++ solver for the current machine.

    This removes any stale precompiled binary, which avoids the common Colab
    "Exec format error" when a macOS executable has accidentally been copied.
    The Bessel cache is not removed.
    """
    global _BUILT_CPP
    if _BUILT_CPP and not force_rebuild:
        return
    if force_rebuild or EXE.exists():
        try:
            EXE.unlink()
        except FileNotFoundError:
            pass
    subprocess.run(["make"], cwd=CPP_DIR, check=True)
    _BUILT_CPP = True


def clone_or_update_repository(
    repo_url: str = "https://github.com/tsmith2/GGI_cosmology_notebooks.git",
    *,
    repo_dir: str | Path = "/content/GGI_cosmology_notebooks",
    branch: str = "main",
) -> Path:
    """Download or update the GitHub repository in a Colab runtime."""
    repo_dir = Path(repo_dir)
    if repo_dir.exists():
        subprocess.run(["git", "-C", str(repo_dir), "fetch", "origin", branch], check=True)
        subprocess.run(["git", "-C", str(repo_dir), "checkout", branch], check=True)
        subprocess.run(["git", "-C", str(repo_dir), "pull", "--ff-only"], check=True)
    else:
        subprocess.run(["git", "clone", "--branch", branch, repo_url, str(repo_dir)], check=True)
    return repo_dir


class TwoFluidCPPServer:
    """Persistent C++ worker for repeated TT spectrum evaluations."""

    def __init__(
        self,
        *,
        ell_max: int = ELL_MAX,
        n_source: int = N_SOURCE,
        k_samples_per_period: float = K_SAMPLES_PER_PERIOD,
        n_threads: int = N_THREADS,
    ) -> None:
        ensure_cpp_executable()
        CACHE_DIR.mkdir(exist_ok=True)
        self.ell_max = int(ell_max)
        self.ell_min = 2
        cmd = [
            str(EXE),
            "--ell-grid",
            "class",
            "--ell-min",
            "2",
            "--ell-max",
            str(self.ell_max),
            "--interpolated-output",
            "--n-source",
            str(int(n_source)),
            "--k-samples-per-period",
            str(float(k_samples_per_period)),
            "--n-threads",
            str(int(n_threads)),
            "--bessel-x-cache",
            str(BESSEL_X_CACHE),
            "--bessel-x-dx",
            "0.10",
            "--bessel-x-max",
            "15000.0",
            "--server",
        ]
        self.process = subprocess.Popen(
            cmd,
            cwd=CPP_DIR,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        ready = self.process.stderr.readline().strip()
        if ready != "SERVER_READY":
            self.close()
            raise RuntimeError(f"C++ worker did not start correctly: {ready}")
        atexit.register(self.close)

    def tt_spectrum(
        self,
        *,
        A_s: float,
        omega_cdm: float,
        omega_b: float,
        H0: float,
        n_s: float,
        Delta_N_eff: float,
    ) -> TTSpectrum:
        """Return the dense TT spectrum in microkelvin squared."""
        if self.process.poll() is not None:
            raise RuntimeError("C++ worker is not running")
        assert self.process.stdin is not None
        assert self.process.stdout is not None

        line = " ".join(
            f"{float(x):.17g}"
            for x in [A_s, omega_cdm, omega_b, H0, n_s, Delta_N_eff]
        )
        self.process.stdin.write(line + "\n")
        self.process.stdin.flush()
        response = self.process.stdout.readline()
        if not response:
            raise RuntimeError("C++ worker closed without a response")
        if response.startswith("ERR "):
            raise RuntimeError(response.strip())
        if not response.startswith("OK "):
            raise RuntimeError(f"bad C++ worker response: {response[:120]}")

        parts = response.split(maxsplit=2)
        n = int(parts[1])
        values = np.fromstring(parts[2], sep=" ")
        if values.size != n:
            raise RuntimeError(f"C++ worker returned {values.size} values, expected {n}")
        ell = np.arange(self.ell_min, self.ell_min + n, dtype=int)
        cl = values * (2.0 * np.pi) / (ell * (ell + 1.0))
        return TTSpectrum(ell=ell, cl=cl, dell_over_2pi=values)

    def close(self) -> None:
        proc = getattr(self, "process", None)
        if proc is None or proc.poll() is not None:
            return
        try:
            if proc.stdin is not None:
                proc.stdin.write("quit\n")
                proc.stdin.flush()
        except Exception:
            pass
        try:
            proc.wait(timeout=2)
        except Exception:
            proc.kill()


def get_server(
    *,
    ell_max: int = ELL_MAX,
    n_source: int = N_SOURCE,
    k_samples_per_period: float = K_SAMPLES_PER_PERIOD,
    n_threads: int = N_THREADS,
) -> TwoFluidCPPServer:
    """Start or reuse a persistent C++ worker."""
    key = (int(ell_max), int(n_source), float(k_samples_per_period), int(n_threads))
    server = _SERVERS.get(key)
    if server is None or server.process.poll() is not None:
        server = TwoFluidCPPServer(
            ell_max=ell_max,
            n_source=n_source,
            k_samples_per_period=k_samples_per_period,
            n_threads=n_threads,
        )
        _SERVERS[key] = server
    return server


def close_servers() -> None:
    """Close any persistent C++ workers."""
    for server in list(_SERVERS.values()):
        server.close()
    _SERVERS.clear()


def load_planck_tt_binning(
    path: str | Path = BINNING_FILE,
    *,
    lmax: int = ELL_MAX,
) -> dict[str, np.ndarray]:
    """Load the compact Planck-like TT binning file included in the package."""
    z = np.load(path)
    blmin = np.array(z["blmin"], dtype=int)
    blmax = np.minimum(np.array(z["blmax"], dtype=int), int(lmax))
    center = np.array(z["center"], dtype=float)
    valid = blmin <= int(lmax)
    return {
        "blmin": blmin[valid],
        "blmax": blmax[valid],
        "center": center[valid],
        "weights_for_dl": np.array(z["weights_for_dl"], dtype=float),
    }


def bin_d_ell(
    d_ell: np.ndarray,
    sigma_d_ell: np.ndarray | None,
    binning: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray | None]:
    """Bin D_ell using the package's Planck-like TT weights."""
    band = []
    sigma = []
    weights = binning["weights_for_dl"]
    for lo, hi in zip(binning["blmin"], binning["blmax"]):
        if weights.size > hi and np.any(weights[lo : hi + 1]):
            w = weights[lo : hi + 1]
        else:
            w = np.ones(hi - lo + 1, dtype=float) / float(hi - lo + 1)
        band.append(float(np.dot(d_ell[lo : hi + 1], w)))
        if sigma_d_ell is not None:
            sigma.append(float(np.sqrt(np.sum((w * sigma_d_ell[lo : hi + 1]) ** 2))))
    return np.array(band), np.array(sigma) if sigma_d_ell is not None else None


def load_mock_data(path: str | Path = MOCK_DATA_FILE) -> dict[str, np.ndarray]:
    """Load the included low-ell plus Planck-binned mock TT data."""
    z = np.load(path)
    return {name: np.array(z[name]) for name in z.files}


def log_prior_lcdm(theta: np.ndarray) -> float:
    """Flat top-hat prior over LCDM_BOUNDS; -inf outside the box."""
    theta = np.asarray(theta, dtype=float)
    if theta.shape != (len(FID_THETA),):
        return -np.inf
    if np.any(theta <= LCDM_BOUNDS[:, 0]) or np.any(theta >= LCDM_BOUNDS[:, 1]):
        return -np.inf
    return 0.0


def _full_tt_spectrum(
    theta: np.ndarray,
    *,
    ell_max: int = ELL_MAX,
) -> TTSpectrum:
    A_s, omega_cdm, omega_b, H0, n_s, Delta_N_eff = map(float, _validate_theta(theta))
    return get_server(ell_max=ell_max).tt_spectrum(
        A_s=A_s,
        omega_cdm=omega_cdm,
        omega_b=omega_b,
        H0=H0,
        n_s=n_s,
        Delta_N_eff=Delta_N_eff,
    )


def tt_spectrum(
    theta: np.ndarray,
    *,
    ell_min: int = ELL_MIN,
    ell_max: int = ELL_MAX,
) -> TTSpectrum:
    """Compute the dense TT spectrum for theta and return the requested ell range."""
    spec = _full_tt_spectrum(theta, ell_max=ell_max)
    mask = (spec.ell >= int(ell_min)) & (spec.ell <= int(ell_max))
    return TTSpectrum(
        ell=spec.ell[mask],
        cl=spec.cl[mask],
        dell_over_2pi=spec.dell_over_2pi[mask],
    )


def _data_vector(
    data: dict[str, np.ndarray],
    *,
    ell_min: int = ELL_MIN,
    ell_max: int = ELL_MAX,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    low_mask = (data["low_ell"] >= ell_min) & (data["low_ell"] <= ell_max)
    bin_mask = (data["bin_center"] >= ell_min) & (data["bin_center"] <= ell_max)
    ell = np.concatenate([data["low_ell"][low_mask], data["bin_center"][bin_mask]])
    obs = np.concatenate(
        [data["low_fake_d_ell"][low_mask], data["fake_d_ell_binned"][bin_mask]]
    )
    sigma = np.concatenate(
        [data["low_sigma_d_ell"][low_mask], data["sigma_d_ell_binned"][bin_mask]]
    )
    return ell, obs, sigma


def model_data_vector(
    theta: np.ndarray,
    data: dict[str, np.ndarray] | None = None,
    *,
    ell_min: int = ELL_MIN,
    ell_max: int = ELL_MAX,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute model values at the same low-ell points and bandpowers as the data."""
    if data is None:
        data = load_mock_data()
    spec = _full_tt_spectrum(theta, ell_max=ell_max)
    full = np.zeros(int(ell_max) + 1, dtype=float)
    in_range = spec.ell <= int(ell_max)
    full[spec.ell[in_range]] = spec.dell_over_2pi[in_range]

    low_mask = (data["low_ell"] >= ell_min) & (data["low_ell"] <= ell_max)
    low_model = full[data["low_ell"][low_mask].astype(int)]

    valid_bins = data["bin_lmin"].astype(int) <= int(ell_max)
    center_all = data["bin_center"].astype(float)[valid_bins]
    binning = {
        "blmin": data["bin_lmin"].astype(int)[valid_bins],
        "blmax": np.minimum(data["bin_lmax"].astype(int)[valid_bins], int(ell_max)),
        "center": center_all,
        "weights_for_dl": load_planck_tt_binning(lmax=ell_max)["weights_for_dl"],
    }
    native_cl, _ = bin_d_ell(full, None, binning)
    center = center_all
    binned_model = native_cl * center * (center + 1.0) / (2.0 * np.pi)
    bin_mask = (center >= ell_min) & (center <= ell_max)

    ell = np.concatenate([data["low_ell"][low_mask], center[bin_mask]])
    model = np.concatenate([low_model, binned_model[bin_mask]])
    return ell, model


def log_likelihood(
    theta: np.ndarray,
    data: dict[str, np.ndarray] | None = None,
    *,
    ell_min: int = ELL_MIN,
    ell_max: int = ELL_MAX,
) -> float:
    """Gaussian log likelihood for the included mock TT data."""
    if data is None:
        data = load_mock_data()
    if not np.isfinite(log_prior_lcdm(theta)):
        return -np.inf
    _ell, obs, sigma = _data_vector(data, ell_min=ell_min, ell_max=ell_max)
    _model_ell, model = model_data_vector(
        theta, data, ell_min=ell_min, ell_max=ell_max
    )
    if np.any(~np.isfinite(model)) or np.any(sigma <= 0.0):
        return -np.inf
    resid = (obs - model) / sigma
    return -0.5 * float(np.sum(resid**2))


def negative_log_likelihood(
    theta: np.ndarray,
    data: dict[str, np.ndarray] | None = None,
    *,
    ell_min: int = ELL_MIN,
    ell_max: int = ELL_MAX,
) -> float:
    """Return -log likelihood for optimizers, scans, or MCMC drivers."""
    ll = log_likelihood(
        theta, data, ell_min=ell_min, ell_max=ell_max
    )
    if not np.isfinite(ll):
        return np.inf
    return -ll


neg_log_likelihood_lcdm = negative_log_likelihood


def plot_power_spectrum(
    theta: np.ndarray = FID_THETA,
    *,
    data: dict[str, np.ndarray] | None = None,
    ell_min: int = ELL_MIN,
    ell_max: int = ELL_MAX,
    output: str | Path | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Plot the model TT power spectrum and, optionally, the mock data."""
    spec = tt_spectrum(theta, ell_min=ell_min, ell_max=ell_max)
    fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=150)
    ax.plot(spec.ell, spec.dell_over_2pi, color="crimson", lw=1.8, label="fiducial model")
    if data is not None:
        low_mask = (data["low_ell"] >= ell_min) & (data["low_ell"] <= ell_max)
        bin_mask = (data["bin_center"] >= ell_min) & (data["bin_center"] <= ell_max)
        ax.errorbar(
            data["low_ell"][low_mask],
            data["low_fake_d_ell"][low_mask],
            yerr=data["low_sigma_d_ell"][low_mask],
            fmt="o",
            ms=2.4,
            lw=0.6,
            color="black",
            ecolor="0.45",
            label="mock low ell",
        )
        ax.errorbar(
            data["bin_center"][bin_mask],
            data["fake_d_ell_binned"][bin_mask],
            yerr=data["sigma_d_ell_binned"][bin_mask],
            fmt="o",
            ms=2.0,
            lw=0.5,
            color="#11a8c8",
            ecolor="#99dbe8",
            alpha=0.9,
            label="mock Planck bins",
        )
    ax.set_xscale("log")
    ax.set_xlim(ell_min, ell_max)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_ylabel(r"$D_\ell^{TT}\;[\mu{\rm K}^2]$")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output)
    return fig, ax


def plot_fiducial_power_spectrum(
    *,
    ell_min: int = ELL_MIN,
    ell_max: int = ELL_MAX,
    include_mock_data: bool = True,
) -> tuple[plt.Figure, plt.Axes]:
    """Convenience plot for the default fiducial parameter vector."""
    data = load_mock_data() if include_mock_data else None
    return plot_power_spectrum(FID_THETA, data=data, ell_min=ell_min, ell_max=ell_max)


def plot_mock_data_with_residuals(
    theta: np.ndarray = FID_THETA,
    *,
    data: dict[str, np.ndarray] | None = None,
    ell_min: int = ELL_MIN,
    ell_max: int = ELL_MAX,
    output: str | Path | None = None,
) -> tuple[plt.Figure, tuple[plt.Axes, plt.Axes]]:
    """Plot mock TT data and residuals normalized by binned cosmic variance."""
    if data is None:
        data = load_mock_data()
    spec = tt_spectrum(theta, ell_min=ell_min, ell_max=ell_max)

    low_mask = (data["low_ell"] >= ell_min) & (data["low_ell"] <= ell_max)
    bin_mask = (data["bin_center"] >= ell_min) & (data["bin_center"] <= ell_max)

    fig, (ax, rx) = plt.subplots(
        2,
        1,
        figsize=(8.2, 6.0),
        dpi=150,
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    ax.plot(spec.ell, spec.dell_over_2pi, color="crimson", lw=1.7, label="model")
    ax.errorbar(
        data["low_ell"][low_mask],
        data["low_fake_d_ell"][low_mask],
        yerr=data["low_sigma_d_ell"][low_mask],
        fmt="o",
        ms=2.3,
        lw=0.6,
        color="black",
        ecolor="0.45",
        label="low ell",
    )
    ax.errorbar(
        data["bin_center"][bin_mask],
        data["fake_d_ell_binned"][bin_mask],
        yerr=data["sigma_d_ell_binned"][bin_mask],
        fmt="o",
        ms=2.0,
        lw=0.5,
        color="#11a8c8",
        ecolor="#99dbe8",
        alpha=0.9,
        label="Planck-like bins",
    )
    ax.set_xscale("log")
    ax.set_xlim(ell_min, ell_max)
    ax.set_ylim(bottom=0.0)
    ax.set_ylabel(r"$D_\ell^{TT}\;[\mu{\rm K}^2]$")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    _ell_model, model = model_data_vector(
        theta, data, ell_min=ell_min, ell_max=ell_max
    )
    _ell_data, obs, sigma = _data_vector(data, ell_min=ell_min, ell_max=ell_max)
    cv_sigma = np.concatenate(
        [
            data["low_cv_sigma_d_ell"][low_mask],
            data["cv_sigma_d_ell_binned"][bin_mask],
        ]
    )
    rx.axhline(0.0, color="crimson", lw=1.0)
    rx.errorbar(
        _ell_data,
        (obs - model) / cv_sigma,
        yerr=sigma / cv_sigma,
        fmt="o",
        ms=2.0,
        lw=0.5,
        color="black",
        ecolor="0.65",
        alpha=0.9,
    )
    rx.set_xlabel(r"Multipole $\ell$")
    rx.set_ylabel(r"$\Delta D_\ell/\sigma_{\rm CV}$")
    rx.set_ylim(-4.0, 4.0)
    rx.grid(alpha=0.25)
    fig.tight_layout()
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output)
    return fig, (ax, rx)


def main() -> None:
    data = load_mock_data()
    nll = negative_log_likelihood(FID_THETA, data)
    FIGURE_DIR.mkdir(exist_ok=True)
    path = FIGURE_DIR / "fiducial_lite_tt_with_mock.png"
    plot_power_spectrum(FID_THETA, data=data, output=path)
    print("theta =", FID_THETA)
    print("parameter names =", PARAM_NAMES)
    print(f"-log likelihood at fiducial theta = {nll:.6g}")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
