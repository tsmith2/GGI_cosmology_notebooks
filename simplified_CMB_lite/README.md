# simplified_CMB_lite

This is a stripped-down RECFAST-enabled CMB temperature solver intended for
fast experiments.  It removes the ISW line-of-sight integral and computes only
the instantaneous last-scattering Sachs-Wolfe plus Doppler TT spectrum.

Baryon density perturbations are kept in the gravitational source, but are not
independently evolved:

```math
\delta_b = {3\over 4}\delta_\gamma .
```

This is approximately true during tight coupling.  The included compact
RECFAST-style module computes the recombination history.  The CMB projection
still assumes instantaneous recombination, but uses RECFAST to compute
`z_star`, the finite-width damping scale, and the Silk damping scale.

The independent mode evolutions are threaded.  By default `--n-threads 0`
uses the number of hardware threads reported by the machine.

Build:

```bash
cd cpp
make
```

Remove a compiled binary:

```bash
make clean
```

This also removes generated Bessel cache files in `cpp/` and `.cache/`.

Clean and rebuild from source:

```bash
make rebuild
```

Run a small smoke spectrum:

```bash
cpp/two_fluid_tt --recfast-output recfast_history.dat \
  --ell-grid sparse --ell-min 2 --ell-max 800 --ell-step 5 \
  --bessel-class-memory --output smoke_tt.dat
```

Run the fast cached path from the package root:

```bash
mkdir -p .cache data
cpp/two_fluid_tt --ell-grid class --ell-min 2 --ell-max 2500 \
  --interpolated-output --n-source 1001 --k-samples-per-period 40 \
  --n-threads 0 \
  --bessel-x-cache .cache/lite_bessel_x_l2500_dx0p1_x15000.bin \
  --bessel-x-dx 0.10 --bessel-x-max 15000 \
  --output data/lite_tt.dat
```

See `explanatory.md` for the physical equations and corresponding C++ line
locations.

## Colab Wrapper And Mock Data

The package includes Planck-like mock TT data generated from the fiducial lite
spectrum:

```text
data/mock_tt_planck_like_lite.npz
data/planck_tt_binning_lmax2500.npz
```

The mock data use individual low-ell points from `ell=2` through `ell=29`, then
Planck-like high-ell bins from `ell=30` to `ell=2500`.  The default effective
experiment is:

```text
beam_fwhm_arcmin = 5.0
noise_uK_arcmin = 45.0
f_sky = 0.70
```

The student-facing wrapper is:

```text
python/simplified_cmb_lite_colab.py
```

The core likelihood call is:

```python
negative_log_likelihood(theta, data=None, ell_min=2, ell_max=2500)
```

where:

```python
theta = [log10(10^9 A_s), omega_cdm, omega_b, h, n_s]
```

The Colab notebook is:

```text
notebooks/simplified_CMB_lite_colab.ipynb
```

It downloads the GitHub repository, compiles the C++ code, plots the fiducial
spectrum, and evaluates the negative log likelihood.  The mock-data residual
plot normalizes `Delta D_ell` by the binned cosmic-variance uncertainty.
