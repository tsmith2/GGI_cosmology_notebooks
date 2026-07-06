# two_fluid_LCDM_colab

This is a slim teaching version of the LCDM two-fluid TT package.  It contains
only the C++ TT solver, a Makefile, one Python wrapper, and the mock
Planck-like TT data needed for a basic likelihood exercise.

The student-facing LCDM parameter vector is

```text
theta = [log10(10^9 A_s), omega_cdm, omega_b, H0, n_s]
```

where `H0` is in km/s/Mpc.  The physical densities `omega_cdm` and `omega_b`
mean `Omega_cdm h^2` and `Omega_b h^2`.

## Compile

From this directory:

```bash
cd cpp
make
cd ..
```

The code also accepts `make NATIVE=1` and `make LTO=1` if you want to try
machine-specific optimization flags.

## Run The Teaching Wrapper

From this directory:

```bash
python python/two_fluid_LCDM_colab.py
```

This builds the C++ executable if needed, computes the fiducial TT spectrum,
loads the included mock data, prints the fiducial `-log likelihood`, and saves

```text
figures/fiducial_lcdm_tt.png
```

The first run may create

```text
cache/two_fluid_bessel_x_class_l2500_dx0p1_x15000_v1.bin
```

which is a reusable Bessel table.

## Notebook Use

In a notebook or Colab cell:

```python
import numpy as np

from python.two_fluid_LCDM_colab import (
    load_mock_data,
    plot_power_spectrum,
    neg_log_likelihood_lcdm,
)

data = load_mock_data()

# Same values as the package's FID_THETA:
# [log10(10^9 A_s), omega_cdm, omega_b, H0, n_s]
theta = np.array([np.log10(2.1), 0.1201, 0.0223, 67.0, 0.965])

print(neg_log_likelihood_lcdm(theta, data))
plot_power_spectrum(theta, data=data)
```

The key functions are:

```text
tt_spectrum(theta)
model_bandpowers(theta)
log_prior_lcdm(theta)
log_likelihood_lcdm(theta, data)
neg_log_likelihood_lcdm(theta, data)
plot_power_spectrum(theta, data=None)
```

No sampler driver is included in this package.
