#!/usr/bin/env python3
"""Helpers for Planck-like mock TT data from a fiducial D_ell spectrum."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


PACKAGE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_BINNING_FILE = PACKAGE_DIR / "data" / "planck_tt_binning_lmax2500.npz"


@dataclass(frozen=True)
class MockTTData:
    ell: np.ndarray
    fid_d_ell: np.ndarray
    noise_d_ell: np.ndarray
    sigma_d_ell: np.ndarray
    cv_sigma_d_ell: np.ndarray
    low_ell: np.ndarray
    low_fid_d_ell: np.ndarray
    low_fake_d_ell: np.ndarray
    low_sigma_d_ell: np.ndarray
    low_cv_sigma_d_ell: np.ndarray
    bin_center: np.ndarray
    bin_lmin: np.ndarray
    bin_lmax: np.ndarray
    fid_native_cl: np.ndarray
    fake_native_cl: np.ndarray
    sigma_native_cl: np.ndarray
    noise_native_cl: np.ndarray
    fid_d_ell_binned: np.ndarray
    fake_d_ell_binned: np.ndarray
    sigma_d_ell_binned: np.ndarray
    cv_sigma_d_ell_binned: np.ndarray
    noise_d_ell_binned: np.ndarray
    f_sky: float
    beam_fwhm_arcmin: float
    noise_uK_arcmin: float
    seed: int


def load_tt_spectrum(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Load ell and D_ell = ell(ell+1) C_ell / 2pi from a text spectrum file."""
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] >= 3:
        return data[:, 0].astype(int), data[:, 2].astype(float)
    if data.shape[1] == 2:
        return data[:, 0].astype(int), data[:, 1].astype(float)
    raise ValueError("TT spectrum file must have at least two columns")


def load_planck_tt_binning(
    path: str | Path = DEFAULT_BINNING_FILE,
    *,
    lmax: int = 2500,
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


def planck_like_noise_d_ell(
    ell: np.ndarray,
    fid_d_ell: np.ndarray,
    *,
    beam_fwhm_arcmin: float = 5.0,
    noise_uK_arcmin: float = 45.0,
    f_sky: float = 0.70,
) -> tuple[np.ndarray, np.ndarray]:
    """Return white-noise D_ell and Gaussian per-ell sigma for a Planck-like TT map."""
    ell = np.asarray(ell, dtype=float)
    fid_d_ell = np.asarray(fid_d_ell, dtype=float)
    sigma_beam = np.deg2rad(beam_fwhm_arcmin / 60.0) / np.sqrt(8.0 * np.log(2.0))
    noise_rad = noise_uK_arcmin * np.pi / (180.0 * 60.0)
    n_ell_cl = noise_rad**2 * np.exp(ell * (ell + 1.0) * sigma_beam**2)
    noise_d_ell = ell * (ell + 1.0) * n_ell_cl / (2.0 * np.pi)
    sigma_d_ell = np.sqrt(2.0 / (f_sky * (2.0 * ell + 1.0))) * (
        fid_d_ell + noise_d_ell
    )
    return noise_d_ell, sigma_d_ell


def bin_d_ell(
    d_ell: np.ndarray,
    sigma_d_ell: np.ndarray | None,
    binning: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray | None]:
    """Bin D_ell using Planck-like weights.

    For native Planck weights, the returned bandpowers are C_ell-like native
    quantities.  Multiply by ell(ell+1)/2pi at the bin center for plotting.
    """
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


def make_mock_tt_data(
    ell: np.ndarray,
    fid_d_ell: np.ndarray,
    *,
    beam_fwhm_arcmin: float = 5.0,
    noise_uK_arcmin: float = 45.0,
    f_sky: float = 0.70,
    seed: int = 20260707,
    lmax: int = 2500,
    binning_file: str | Path = DEFAULT_BINNING_FILE,
) -> MockTTData:
    """Generate low-ell plus Planck-binned high-ell mock TT data."""
    ell = np.asarray(ell, dtype=int)
    fid_d_ell = np.asarray(fid_d_ell, dtype=float)
    if ell.shape != fid_d_ell.shape:
        raise ValueError("ell and fid_d_ell must have the same shape")
    if np.min(ell) < 2:
        raise ValueError("ell values must start at ell >= 2")

    lmax = min(int(lmax), int(np.max(ell)))
    keep = ell <= lmax
    ell = ell[keep]
    fid_d_ell = fid_d_ell[keep]
    binning = load_planck_tt_binning(binning_file, lmax=lmax)

    full_ell = np.arange(lmax + 1, dtype=float)
    full_d_ell = np.zeros(lmax + 1, dtype=float)
    full_d_ell[ell] = fid_d_ell

    noise_d_ell, sigma_d_ell = planck_like_noise_d_ell(
        ell.astype(float),
        fid_d_ell,
        beam_fwhm_arcmin=beam_fwhm_arcmin,
        noise_uK_arcmin=noise_uK_arcmin,
        f_sky=f_sky,
    )
    cv_sigma_d_ell = np.sqrt(2.0 / (f_sky * (2.0 * ell.astype(float) + 1.0))) * fid_d_ell
    full_noise_d_ell = np.zeros_like(full_ell)
    full_sigma_d_ell = np.zeros_like(full_ell)
    full_cv_sigma_d_ell = np.zeros_like(full_ell)
    full_noise_d_ell[ell] = noise_d_ell
    full_sigma_d_ell[ell] = sigma_d_ell
    full_cv_sigma_d_ell[ell] = cv_sigma_d_ell

    fid_native, sigma_native = bin_d_ell(full_d_ell, full_sigma_d_ell, binning)
    _fid_native_cv, cv_sigma_native = bin_d_ell(full_d_ell, full_cv_sigma_d_ell, binning)
    noise_native, _ = bin_d_ell(full_noise_d_ell, None, binning)
    rng = np.random.default_rng(seed)
    fake_native = fid_native + rng.normal(0.0, sigma_native)

    first_bin = int(np.min(binning["blmin"]))
    low_mask = ell < first_bin
    low_ell = ell[low_mask]
    low_fid = fid_d_ell[low_mask]
    low_sigma = sigma_d_ell[low_mask]
    low_cv_sigma = cv_sigma_d_ell[low_mask]
    low_fake = low_fid + rng.normal(0.0, low_sigma)

    center = binning["center"]
    cl_to_dl_center = center * (center + 1.0) / (2.0 * np.pi)
    fid_plot = fid_native * cl_to_dl_center
    fake_plot = fake_native * cl_to_dl_center
    sigma_plot = sigma_native * cl_to_dl_center
    cv_sigma_plot = cv_sigma_native * cl_to_dl_center
    noise_plot = noise_native * cl_to_dl_center

    return MockTTData(
        ell=ell,
        fid_d_ell=fid_d_ell,
        noise_d_ell=noise_d_ell,
        sigma_d_ell=sigma_d_ell,
        cv_sigma_d_ell=cv_sigma_d_ell,
        low_ell=low_ell,
        low_fid_d_ell=low_fid,
        low_fake_d_ell=low_fake,
        low_sigma_d_ell=low_sigma,
        low_cv_sigma_d_ell=low_cv_sigma,
        bin_center=center,
        bin_lmin=binning["blmin"],
        bin_lmax=binning["blmax"],
        fid_native_cl=fid_native,
        fake_native_cl=fake_native,
        sigma_native_cl=sigma_native,
        noise_native_cl=noise_native,
        fid_d_ell_binned=fid_plot,
        fake_d_ell_binned=fake_plot,
        sigma_d_ell_binned=sigma_plot,
        cv_sigma_d_ell_binned=cv_sigma_plot,
        noise_d_ell_binned=noise_plot,
        f_sky=float(f_sky),
        beam_fwhm_arcmin=float(beam_fwhm_arcmin),
        noise_uK_arcmin=float(noise_uK_arcmin),
        seed=int(seed),
    )


def save_mock_tt_data(mock: MockTTData, path: str | Path) -> None:
    """Save a mock TT data set to an npz file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **mock.__dict__)
