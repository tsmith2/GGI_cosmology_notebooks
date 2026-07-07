#!/usr/bin/env python3
"""Teaching wrapper for the lite RECFAST/theta two-fluid TT solver.

Parameter vector:

    theta = [log10(10^9 A_s), omega_cdm, omega_b, h, n_s]

where omega_cdm and omega_b are physical densities, Omega_i h^2.
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
MPLCONFIG_DIR = PACKAGE_DIR / ".mplconfig"
XDG_CACHE_DIR = PACKAGE_DIR / ".cache"
MPLCONFIG_DIR.mkdir(exist_ok=True)
XDG_CACHE_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(XDG_CACHE_DIR))

import matplotlib.pyplot as plt

ELL_MIN = 2
ELL_MAX = 2500
N_SOURCE = 1001
K_SAMPLES_PER_PERIOD = 40.0
N_THREADS = 0
BESSEL_X_CACHE = PACKAGE_DIR / ".cache" / "lite_bessel_x_l2500_dx0p1_x15000.bin"

FID_THETA = np.array([np.log10(2.1), 0.1201, 0.0223, 0.67, 0.965], dtype=float)
_BUILT_CPP = False


@dataclass(frozen=True)
class TTSpectrum:
    ell: np.ndarray
    cl: np.ndarray
    dell_over_2pi: np.ndarray


def As_from_log10_1e9_As(log10_1e9_As: float) -> float:
    """Convert log10(10^9 A_s) into A_s."""
    return 10.0 ** float(log10_1e9_As) * 1.0e-9


def ensure_cpp_executable() -> None:
    """Build the C++ solver for the current machine.

    The repository should not rely on a precompiled executable.  In particular,
    a Mac binary copied into GitHub will fail in Colab with "Exec format error",
    so we remove only the binary before compiling.  The fixed-x Bessel cache is
    intentionally left in place if it already exists.
    """
    global _BUILT_CPP
    if _BUILT_CPP:
        return
    try:
        EXE.unlink()
    except FileNotFoundError:
        pass
    subprocess.run(["make"], cwd=CPP_DIR, check=True)
    _BUILT_CPP = True


class TwoFluidCPPServer:
    """Persistent C++ worker for the lite primary TT spectrum."""

    def __init__(
        self,
        *,
        ell_min: int = ELL_MIN,
        ell_max: int = ELL_MAX,
        n_source: int = N_SOURCE,
        k_samples_per_period: float = K_SAMPLES_PER_PERIOD,
        n_threads: int = N_THREADS,
        bessel_x_cache: Path | str = BESSEL_X_CACHE,
    ) -> None:
        ensure_cpp_executable()
        self.ell_min = int(ell_min)
        self.ell_max = int(ell_max)
        cmd = [
            str(EXE),
            "--ell-grid",
            "class",
            "--ell-min",
            str(self.ell_min),
            "--ell-max",
            str(self.ell_max),
            "--interpolated-output",
            "--bessel-x-cache",
            str(bessel_x_cache),
            "--n-source",
            str(int(n_source)),
            "--k-samples-per-period",
            str(float(k_samples_per_period)),
            "--n-threads",
            str(int(n_threads)),
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
        h: float,
        n_s: float,
    ) -> TTSpectrum:
        """Return the TT spectrum in microkelvin squared."""
        if self.process.poll() is not None:
            raise RuntimeError("C++ worker is not running")
        assert self.process.stdin is not None
        assert self.process.stdout is not None

        line = " ".join(
            f"{float(x):.17g}" for x in [A_s, omega_cdm, omega_b, h, n_s]
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
        return TTSpectrum(
            ell=ell,
            cl=cl,
            dell_over_2pi=values,
        )

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


def tt_spectrum(theta: np.ndarray) -> TTSpectrum:
    """Compute one TT spectrum from theta."""
    log_as, omega_cdm, omega_b, h, n_s = map(float, theta)
    server = TwoFluidCPPServer()
    try:
        return server.tt_spectrum(
            A_s=As_from_log10_1e9_As(log_as),
            omega_cdm=omega_cdm,
            omega_b=omega_b,
            h=h,
            n_s=n_s,
        )
    finally:
        server.close()


def plot_tt_spectrum(
    spec: TTSpectrum,
    *,
    xscale: str = "log",
) -> tuple[plt.Figure, plt.Axes]:
    """Plot D_ell = ell(ell+1) C_ell / 2pi for the lite TT spectrum."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(spec.ell, spec.dell_over_2pi, lw=2.0, color="black")
    ax.set_xlim(2, 2500)
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_ylabel(r"$\ell(\ell+1)C_\ell^{TT}/2\pi\;[\mu{\rm K}^2]$")
    ax.set_xscale(xscale)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    return fig, ax
