# Key Equations In The Two-Fluid LCDM Code

This file maps the main physics equations in `cpp/two_fluid_tt.cpp` to the
lines where they are implemented.  The code is a simplified instantaneous
recombination, tight-coupling, no-neutrino two-fluid model for the CMB TT power
spectrum.  It is meant for exploration and teaching, not as a replacement for
CLASS or CAMB.

The independent time variable in the perturbation equations is

```math
x = \eta / \tau_r ,
```

and the wavenumber variable is

```math
\kappa = k \tau_r .
```

## Input Parameters And Background

| Idea | Equation / convention | Code location |
|---|---|---|
| Parameter vector | `theta = [log10(10^9 A_s), omega_cdm, omega_b, H0, n_s]`, with `H0` in km/s/Mpc | `Params` defaults at `cpp/two_fluid_tt.cpp:25-31`; command-line parsing at `cpp/two_fluid_tt.cpp:1121-1126`; server input at `cpp/two_fluid_tt.cpp:1246-1248` |
| Convert `H0` to `h` | `h = H0 / 100` | `cpp/two_fluid_tt.cpp:55-57` |
| Convert physical densities to fractional densities | `Omega_cdm = omega_cdm / h^2`, `Omega_b = omega_b / h^2`, `Omega_gamma = 2.47e-5 / h^2` | `cpp/two_fluid_tt.cpp:55-59` |
| Strict flatness | `Omega_Lambda = 1 - Omega_cdm - Omega_b - Omega_rad` | `cpp/two_fluid_tt.cpp:60-62` |
| Conformal time integral | `eta(a) = int_0^a da / [H0 sqrt(Omega_m a + Omega_r + Omega_Lambda a^4)]` | integrand at `cpp/two_fluid_tt.cpp:77-83`; Simpson rule at `cpp/two_fluid_tt.cpp:84-88` |
| Recombination redshift fit | `z_rec = 1000 omega_b^{-0.027/(1+0.11 ln omega_b)}` | `cpp/two_fluid_tt.cpp:327-329` |
| Equality and recombination scales | `a_rec`, `a_eq`, `eta_0`, `eta_rec`, `tau_r`, `x_rec` | `cpp/two_fluid_tt.cpp:143-158` |
| Analytic scale-factor variable | `y(x) = alpha^2 x^2 + 2 alpha x` | `cpp/two_fluid_tt.cpp:331` |

## Perturbation State

The ODE state vector is

```text
u = (delta_gamma, delta_c, v_gamma, v_c, Phi)
```

and is declared as a five-component array at `cpp/two_fluid_tt.cpp:70`.

## Dynamical Equations

The perturbation equations are implemented in `rhs(...)`, beginning at
`cpp/two_fluid_tt.cpp:399`.

| Quantity | Equation / meaning | Code location |
|---|---|---|
| Dimensionless conformal Hubble factor | `eta = 2 alpha (alpha x + 1) / y` in the code's dimensionless variables | `cpp/two_fluid_tt.cpp:404-405` |
| Baryon and CDM loading factors | `y_b = (Omega_b/Omega_m) y`, `y_c = (Omega_cdm/Omega_m) y` | `cpp/two_fluid_tt.cpp:406-407` |
| Newtonian potential evolution | `Phi' = -eta Phi + [3 eta^2 /(2(1+y) kappa)] [(4/3 + y - y_c) v_gamma + y_c v_c]` | `cpp/two_fluid_tt.cpp:408-412` |
| Photon density equation | `delta_gamma' = -(4/3) kappa v_gamma + 4 Phi'` | `cpp/two_fluid_tt.cpp:414` |
| CDM density equation | `delta_c' = -kappa v_c + 3 Phi'` | `cpp/two_fluid_tt.cpp:415` |
| Tightly-coupled photon-baryon velocity | `v_gamma' = [-eta y_b v_gamma + kappa delta_gamma/3]/(4/3+y_b) + kappa Phi` | `cpp/two_fluid_tt.cpp:416-418` |
| CDM velocity | `v_c' = -eta v_c + kappa Phi` | `cpp/two_fluid_tt.cpp:419` |
| Potential state derivative | stored as `u[4]' = Phi'` | `cpp/two_fluid_tt.cpp:420` |

## Initial Conditions

The code uses the conformal-Newtonian growing-mode initial conditions based on
Ma & Bertschinger Eq. 98, with the CLASS-like convention `C = 1/2`.

| Quantity | Equation / convention | Code location |
|---|---|---|
| Initial time | `x_i = min(1e-4/kappa, x_eq/10000)` | `cpp/two_fluid_tt.cpp:424-425` |
| Curvature convention | `C = 1/2` | `cpp/two_fluid_tt.cpp:426-429` |
| Initial potential | `Phi_i = 20 C / 15` | `cpp/two_fluid_tt.cpp:428-429` |
| Photon perturbation | `delta_gamma,i = -2 Phi_i (1 + 3 y_i/16)` | `cpp/two_fluid_tt.cpp:430-432` |
| CDM perturbation | `delta_c,i = 3 delta_gamma,i / 4` | `cpp/two_fluid_tt.cpp:433` |
| Initial velocity | shared photon/CDM velocity `v_gamma,i = v_c,i` from the growing-mode expansion | `cpp/two_fluid_tt.cpp:434-439` |

## ODE Solver

| Idea | Equation / convention | Code location |
|---|---|---|
| Adaptive RK45 step | Dormand-Prince 5(4) coefficients | coefficients at `cpp/two_fluid_tt.cpp:442-469`; stage evaluations at `cpp/two_fluid_tt.cpp:480-505` |
| Error control | `err_norm = max_i |err_i|/(atol + rtol max(|u_i|, |trial_i|))` | `cpp/two_fluid_tt.cpp:520-523` |
| Step-size update | `h -> h * clamp(0.9 err_norm^{-0.2}, 0.2, 5)` | `cpp/two_fluid_tt.cpp:529-531` |
| Integrate each `kappa` to recombination | evolve from `x_i` to `x_rec` | `cpp/two_fluid_tt.cpp:508-535` |

## Temperature Sources At Recombination

The instantaneous-recombination source is built after integrating each mode to
`x_rec`.

| Source | Equation / meaning | Code location |
|---|---|---|
| Sachs-Wolfe source | `S_SW = Phi + delta_gamma/4` | `cpp/two_fluid_tt.cpp:623-637` |
| Doppler source | `S_Doppler = v_gamma` | `cpp/two_fluid_tt.cpp:623-637` |
| Potential at recombination | `Phi_rec = Phi(x_rec)` | `cpp/two_fluid_tt.cpp:638` |
| Radiation fraction proxy | used by the approximate early/late ISW split | `cpp/two_fluid_tt.cpp:640-647` |
| Interpolation of sources in `kappa` | cubic splines for SW, Doppler, `Phi_rec`, and radiation fraction | `cpp/two_fluid_tt.cpp:167-173` |

## Projection To The Angular Power Spectrum

For each multipole, the code computes a line-of-sight-projected source and
integrates over `kappa`.

| Idea | Equation / convention | Code location |
|---|---|---|
| Distance to last scattering | `chi = (eta_0 - eta_*) / tau_r` | `cpp/two_fluid_tt.cpp:183` |
| Bessel argument | `x = kappa chi` | `cpp/two_fluid_tt.cpp:202-205` |
| Bessel derivative | `j_l'(x) = -j_l(x)/(2x) + [j_{l-1}(x)-j_{l+1}(x)]/2` | `cpp/two_fluid_tt.cpp:281-286` |
| Primary source | `Theta_l(kappa) = S_SW j_l(kappa chi) + v_gamma j_l'(kappa chi)` | `cpp/two_fluid_tt.cpp:281-287` |
| Power-spectrum integral | `C_l = 4 pi A_s T_CMB^2 tilt int d kappa [source^2/kappa]` | amplitude at `cpp/two_fluid_tt.cpp:194-197`; integral at `cpp/two_fluid_tt.cpp:279-317` |
| Tilt around pivot | `(l/750)^(n_s - 0.965)` | `cpp/two_fluid_tt.cpp:315` |
| Output `D_l` | `D_l = l(l+1) C_l / (2 pi)` | `cpp/two_fluid_tt.cpp:1189-1204`; file output at `cpp/two_fluid_tt.cpp:1207-1237` |

## Damping Terms

| Term | Equation / convention | Code location |
|---|---|---|
| Sound speed at recombination | `c_s = [3(1 + 3 y_b/4)]^{-1/2}` | `cpp/two_fluid_tt.cpp:339-342` |
| Silk damping scale proxy | `x_s = 0.6 Omega_m^{1/4} Omega_b^{-1/2} a_rec^{3/4} h^{-1/2}` | `cpp/two_fluid_tt.cpp:157-158` |
| Finite visibility width | `sigma_x = (sigma_eta/eta_rec) x_rec`; damping contribution `c_s^2 sigma_x^2` | `cpp/two_fluid_tt.cpp:198-200` |
| Total primary damping | `exp[-kappa^2 (2 x_s^2 + c_s^2 sigma_x^2)]` | `cpp/two_fluid_tt.cpp:204-207` |

## ISW Approximations

The default mode is `early-late`.

| Mode | Equation / convention | Code location |
|---|---|---|
| Seljak-style analytic `Delta Phi` | `Delta Phi = [2 - 8/y_rec + 16 x_rec/y_rec^3]/(10 y_rec)` | `cpp/two_fluid_tt.cpp:333-337`; inserted at `cpp/two_fluid_tt.cpp:188` and `cpp/two_fluid_tt.cpp:208` |
| Early ISW approximation | `ISW_e ~ -2 Phi_* f_r [j_l + (kappa x_rec/3) j_l']` | `cpp/two_fluid_tt.cpp:288-296` |
| Late ISW growth approximation | integrates `d ln D/d ln a = Omega_m(a)^0.55` and compresses the source to an effective distance | `cpp/two_fluid_tt.cpp:355-397` |
| Late ISW contribution | `ISW_late ~ 2 Phi_* f_m Delta(D/a) j_l(kappa chi_eff)` for low multipoles | `cpp/two_fluid_tt.cpp:297-301`; late Bessel setup at `cpp/two_fluid_tt.cpp:244-267` |
| Full LOS ISW option | directly integrates `2 Phi'(eta) j_l[k(eta_0-eta)]` from recombination to today | `cpp/two_fluid_tt.cpp:538-620` |

## Bessel Functions And Sampling

| Idea | Equation / convention | Code location |
|---|---|---|
| CLASS-like multipole sampling | log steps at low `l`, linear step 40 at high `l` | `cpp/two_fluid_tt.cpp:680-712` |
| Required Bessel orders | need `l-1`, `l`, `l+1` because of `j_l'` | `cpp/two_fluid_tt.cpp:714-722` |
| Spherical Bessel recurrence | downward recurrence normalized to exact `j_0` or `j_1` | `cpp/two_fluid_tt.cpp:724-766` |
| Fixed-`x` Bessel cache | load once into memory and interpolate to requested `kappa chi` values | read/interpolate at `cpp/two_fluid_tt.cpp:994-1065`; write at `cpp/two_fluid_tt.cpp:1067-1110` |

## Python Wrapper Connection

The wrapper `python/two_fluid_LCDM_colab.py` talks to the C++ code in server
mode.  The C++ server reads one parameter line

```text
A_s omega_cdm omega_b H0 n_s
```

and returns the dense `D_l` array.  This is handled at
`cpp/two_fluid_tt.cpp:1239-1267`.
