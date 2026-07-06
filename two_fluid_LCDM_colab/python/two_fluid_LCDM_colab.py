#!/usr/bin/env python3
"""Small teaching wrapper for the two-fluid LCDM TT code.

The functions here are meant to be copied into, or imported from, a notebook.
They keep the same likelihood structure as the MCMC driver from the larger
package, but remove the sampler-specific machinery.

Parameter vector:

    theta = [log10(10^9 A_s), omega_cdm, omega_b, H0, n_s]

where H0 is in km/s/Mpc, and omega_cdm and omega_b are physical densities
Omega_cdm h^2 and Omega_b h^2.
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
CACHE_DIR = PACKAGE_DIR / "cache"
FIGURE_DIR = PACKAGE_DIR / "figures"
MOCK_DATA_FILE = PACKAGE_DIR / "data" / "mock_tt_planck_like.npz"
BINNING_FILE = PACKAGE_DIR / "data" / "planck_tt_binning_lmax2500.npz"
BESSEL_X_CACHE = CACHE_DIR / "two_fluid_bessel_x_class_l2500_dx0p1_x15000_v1.bin"

MPLCONFIG_DIR = PACKAGE_DIR / ".mplconfig"
XDG_CACHE_DIR = PACKAGE_DIR / ".cache"
MPLCONFIG_DIR.mkdir(exist_ok=True)
XDG_CACHE_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(XDG_CACHE_DIR))

import matplotlib.pyplot as plt

PARAM_NAMES = [
    r"$\log_{10}(10^9 A_s)$",
    r"$\omega_{\rm cdm}$",
    r"$\omega_b$",
    r"$H_0$",
    r"$n_s$",
]

FID_THETA = np.array([np.log10(2.1), 0.1201, 0.0223, 67.0, 0.965], dtype=float)
LCDM_BOUNDS = np.array(
    [
        [-0.2, 0.8],
        [0.09, 0.16],
        [0.017, 0.028],
        [55.0, 80.0],
        [0.85, 1.08],
    ],
    dtype=float,
)

ELL_MIN = 2
ELL_MAX = 2500
N_SOURCE = 1001
SERVER: TwoFluidCPPServer | None = None


@dataclass(frozen=True)
class TTSpectrum:
    ell: np.ndarray
    cl: np.ndarray
    dell_over_2pi: np.ndarray


def As_from_log10_1e9_As(log10_1e9_As: float) -> float:
    """Convert log10(10^9 A_s) into A_s."""
    return 10.0 ** float(log10_1e9_As) * 1.0e-9


def ensure_cpp_executable() -> None:
    """Build the C++ solver if it is not already present."""
    if EXE.exists():
        return
    subprocess.run(["make"], cwd=CPP_DIR, check=True)


class TwoFluidCPPServer:
    """Persistent C++ worker for repeated likelihood evaluations."""

    def __init__(
        self,
        *,
        ell_min: int = ELL_MIN,
        ell_max: int = ELL_MAX,
        n_source: int = N_SOURCE,
    ) -> None:
        ensure_cpp_executable()
        CACHE_DIR.mkdir(exist_ok=True)
        self.ell_min = int(ell_min)
        self.ell_max = int(ell_max)
        cmd = [
            str(EXE),
            "--n-source",
            str(int(n_source)),
            "--ell-grid",
            "class",
            "--ell-min",
            str(self.ell_min),
            "--ell-max",
            str(self.ell_max),
            "--interpolated-output",
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
    ) -> TTSpectrum:
        """Return D_ell = ell(ell+1) C_ell / 2pi from the C++ solver."""
        if self.process.poll() is not None:
            raise RuntimeError("C++ worker is not running")
        assert self.process.stdin is not None
        assert self.process.stdout is not None
        line = " ".join(
            f"{float(x):.17g}" for x in [A_s, omega_cdm, omega_b, H0, n_s]
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


def get_server() -> TwoFluidCPPServer:
    """Start or reuse the persistent C++ worker."""
    global SERVER
    if SERVER is None or SERVER.process.poll() is not None:
        SERVER = TwoFluidCPPServer()
    return SERVER


def load_planck_tt_binning(path: str | Path = BINNING_FILE, *, lmax: int = ELL_MAX) -> dict[str, np.ndarray]:
    """Load the compact Planck-like TT binning file included in the package."""
    z = np.load(path)
    blmin = np.array(z["blmin"], dtype=int)
    blmax = np.minimum(np.array(z["blmax"], dtype=int), lmax)
    center = np.array(z["center"], dtype=float)
    valid = blmin <= lmax
    return {
        "blmin": blmin[valid],
        "blmax": blmax[valid],
        "center": center[valid],
        "weights_for_dl": np.array(z["weights_for_dl"], dtype=float),
    }


BINNING = load_planck_tt_binning()


def bin_d_ell(
    d_ell: np.ndarray,
    sigma_d_ell: np.ndarray | None,
    binning: dict[str, np.ndarray] = BINNING,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Bin a D_ell array using the package's Planck-like TT bins."""
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
    """Load the included Planck-like mock data with convenient likelihood names."""
    z = np.load(path)
    data = {name: np.array(z[name]) for name in z.files}
    data["d_obs"] = data["fake_d_ell_binned"]
    data["d_sigma"] = data["sigma_d_ell_binned"]
    data["center"] = data["bin_center"]
    return data


def log_prior_lcdm(theta: np.ndarray) -> float:
    """Flat top-hat prior over LCDM_BOUNDS; -inf outside the box."""
    theta = np.asarray(theta, dtype=float)
    if theta.shape != (len(FID_THETA),):
        return -np.inf
    if np.any(theta <= LCDM_BOUNDS[:, 0]) or np.any(theta >= LCDM_BOUNDS[:, 1]):
        return -np.inf
    return 0.0


def tt_spectrum(theta: np.ndarray) -> TTSpectrum:
    """Compute the unbinned TT spectrum for theta."""
    log_as, omega_cdm, omega_b, H0, n_s = map(float, theta)
    return get_server().tt_spectrum(
        A_s=As_from_log10_1e9_As(log_as),
        omega_cdm=omega_cdm,
        omega_b=omega_b,
        H0=H0,
        n_s=n_s,
    )


def model_bandpowers(theta: np.ndarray) -> np.ndarray:
    """Run the C++ solver at theta and bin its TT spectrum into data bandpowers."""
    spec = tt_spectrum(theta)
    full = np.zeros(ELL_MAX + 1)
    full[spec.ell.astype(int)] = spec.dell_over_2pi
    native_cl, _ = bin_d_ell(full, None, BINNING)
    c = BINNING["center"]
    return native_cl * c * (c + 1.0) / (2.0 * np.pi)


def log_likelihood_lcdm(theta: np.ndarray, data: dict[str, np.ndarray]) -> float:
    """Gaussian bandpower log likelihood between the model and mock data."""
    model = model_bandpowers(theta)
    if np.any(~np.isfinite(model)) or np.any(model <= 0.0):
        return -np.inf
    resid = (data["d_obs"] - model) / data["d_sigma"]
    return -0.5 * float(np.sum(resid**2))


def neg_log_likelihood_lcdm(theta: np.ndarray, data: dict[str, np.ndarray]) -> float:
    """Return -log likelihood, useful for optimizers and simple grid searches."""
    ll = log_likelihood_lcdm(theta, data)
    if not np.isfinite(ll):
        return np.inf
    return -ll


def plot_power_spectrum(
    theta: np.ndarray = FID_THETA,
    *,
    data: dict[str, np.ndarray] | None = None,
    output: str | Path | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Plot the model TT power spectrum, optionally over the included mock data."""
    spec = tt_spectrum(theta)
    fig, ax = plt.subplots(figsize=(7.0, 4.4), dpi=160)
    ax.plot(spec.ell, spec.dell_over_2pi, color="crimson", lw=1.6, label="two-fluid LCDM")
    if data is not None:
        ax.errorbar(
            data["center"],
            data["d_obs"],
            yerr=data["d_sigma"],
            fmt="o",
            ms=2.4,
            lw=0.6,
            color="black",
            ecolor="0.45",
            alpha=0.85,
            label="mock bandpowers",
        )
    ax.set_xlim(2, ELL_MAX)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel(r"$\mathcal{D}_\ell^{TT}\,[\mu{\rm K}^2]$")
    ax.legend(frameon=False)
    fig.tight_layout()
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output)
    return fig, ax


# Example "cell" for a notebook:
#
# data = load_mock_data()
# theta = FID_THETA.copy()
# print("negative log likelihood =", neg_log_likelihood_lcdm(theta, data))
# plot_power_spectrum(theta, data=data)


def main() -> None:
    data = load_mock_data()
    nll = neg_log_likelihood_lcdm(FID_THETA, data)
    FIGURE_DIR.mkdir(exist_ok=True)
    plot_path = FIGURE_DIR / "fiducial_lcdm_tt.png"
    plot_power_spectrum(FID_THETA, data=data, output=plot_path)
    print("theta =", FID_THETA)
    print("parameter names =", PARAM_NAMES)
    print(f"-log likelihood at fiducial theta = {nll:.6g}")
    print(f"wrote {plot_path}")


if __name__ == "__main__":
    main()
