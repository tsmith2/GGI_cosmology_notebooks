# Simplified CMB Two-Fluid RECFAST Version

This package evolves a simplified flat LCDM photon-baryon-CDM model in
Newtonian gauge.  The independent time variable is conformal time `eta`, in
Mpc, and the wavenumber is the physical `k`, in Mpc^-1.

The perturbation state is

```text
u = (delta_gamma, delta_b, delta_c, theta_gammab, theta_c, Phi).
```

In this simplified version, `delta_b` is not an independent dynamical
variable.  It is set algebraically from the photon perturbation:

```math
\delta_b = {3\over 4}\delta_\gamma .
```

This is approximately true during tight coupling.  The perturbation state
actually evolved by the code is therefore

```text
u = (delta_gamma, delta_c, theta_gammab, theta_c, Phi).
```

There are no velocity variables `v_gamma` or `v_c` in the code.  The code uses
velocity divergences `theta_gammab` and `theta_c`.

Code location:

- `cpp/two_fluid_tt.cpp:239-240` defines the five-component state.
- `cpp/two_fluid_tt.cpp:748-755` unpacks the state and sets
  `delta_b = 3 delta_gamma / 4` inside the perturbation
  right-hand side.

## Background

The density parameters are built from the supplied physical densities
`omega_cdm = Omega_c h^2`, `omega_b = Omega_b h^2`, and `h`.  The model is
strictly flat:

```math
\Omega_\Lambda = 1 - \Omega_c - \Omega_b - \Omega_\gamma .
```

The photon density is fixed by the CMB temperature convention used in the code,

```math
\Omega_\gamma h^2 = 2.47\times 10^{-5}.
```

Code location:

- `cpp/two_fluid_tt.cpp:65-78` converts physical densities to fractional
  densities and imposes flatness.

The background conformal time is computed from the Friedmann equation,

```math
\left({a'\over a}\right)^2
=
a^2 H_0^2
\left(
{\Omega_\gamma\over a^4}
+ {\Omega_b\over a^3}
+ {\Omega_c\over a^3}
+ \Omega_\Lambda
\right),
```

or equivalently

```math
{d\eta\over d\ln a}
=
{1\over a H(a)}.
```

Code location:

- `cpp/two_fluid_tt.cpp:368-429` integrates the background ODE for
  `eta(a)`.
- `cpp/two_fluid_tt.cpp:476-488` initializes `eta_today`, `eta_rec`,
  `eta_eq`, and the recombination distance.

The compact RECFAST-style module computes `x_e(z)` and `T_m(z)`.  The TT
projection still treats recombination as instantaneous, but the RECFAST history
is used to build the visibility function, choose `z_star`, estimate the
finite-width damping scale, and compute the Silk damping scale.

Code location:

- `cpp/two_fluid_tt.cpp:81-237` contains the compact RECFAST module.
- `cpp/two_fluid_tt.cpp:476-488` computes the RECFAST history and derives
  `z_star`, `sigma_eta`, and the Silk scale.
- `cpp/two_fluid_tt.cpp:496-505` writes the RECFAST history to an ASCII file.
- `cpp/two_fluid_tt.cpp:659-746` builds the visibility function and evaluates
  the diffusion damping integral.

## Perturbations

The potential is evolved by solving the Newtonian-gauge Einstein constraint for
`Phi'`:

```math
k^2\Phi
+ 3{\cal H}\left(\Phi' + {\cal H}\Phi\right)
=
-{3\over 2}{\cal H}^2 \delta_{\rm tot},
```

where

```math
{\cal H} = {a'\over a}
```

and

```math
\delta_{\rm tot}
=
{\rho_\gamma\delta_\gamma
+ \rho_b\delta_b
+ \rho_c\delta_c
\over
\rho_\gamma+\rho_b+\rho_c+\rho_\Lambda}.
```

Here the baryon perturbation is included in the gravitational source, but is
not separately evolved:

```math
\delta_b = {3\over 4}\delta_\gamma .
```

The photon-baryon inertia ratio is

```math
R = {4\rho_\gamma\over 3\rho_b}.
```

The perturbation equations are

```math
\delta_\gamma'
=
-{4\over 3}\theta_{\gamma b}
+ 4\Phi',
```

```math
\delta_c'
=
-\theta_c
+ 3\Phi',
```

```math
\theta_{\gamma b}'
=
-{{\cal H}\over 1+R}\theta_{\gamma b}
+ k^2\Phi
+ {R\over 4(1+R)} k^2\delta_\gamma,
```

```math
\theta_c'
=
-{\cal H}\theta_c
+ k^2\Phi.
```

Code location:

- `cpp/two_fluid_tt.cpp:756-764` computes `a`, densities, `H(a)`, and
  `calH`.
- `cpp/two_fluid_tt.cpp:754-771` sets `delta_b = 3 delta_gamma / 4`, computes
  `delta_tot`, and solves the
  constraint for `Phi'`.
- `cpp/two_fluid_tt.cpp:772-778` evaluates the four fluid equations and
  stores `Phi'`.

## Initial Conditions

For each wavenumber, the integration starts at

```math
\eta_i = \min(10^{-3}/k,\;10^{-3}\eta_{\rm eq}).
```

The initial conditions are

```math
\Phi_i = -{2\over 3},
```

```math
\delta_{\gamma,i} = -2\Phi_i = {4\over 3},
```

```math
\delta_{c,i}
=
{3\over 4}\delta_{\gamma,i}
= 1.
```

The baryon perturbation is not stored in the state; whenever it is needed in
the potential source it is set to

```math
\delta_{b,i} = {3\over 4}\delta_{\gamma,i}=1.
```

```math
\theta_{\gamma b,i}
=
\theta_{c,i}
=
{1\over 2} k^2\eta_i\Phi_i
=
-{1\over 3}k^2\eta_i.
```

The choice `Phi_i = -2/3` is set so that, on superhorizon scales during
radiation domination, `Phi_i = \mathcal{R}`, the primordial curvature
perturbation, in the convention used by this pedagogical code.

Code location:

- `cpp/two_fluid_tt.cpp:782-788` sets the initial time and initial state.

## Temperature Sources

At instantaneous recombination the Sachs-Wolfe source is

```math
S_{\rm SW}(k)
=
\Phi(k,\eta_\ast)
+ {1\over 4}\delta_\gamma(k,\eta_\ast).
```

The Doppler source is written using the photon-baryon velocity divergence:

```math
S_{\rm Doppler}(k)
=
-{\theta_{\gamma b}(k,\eta_\ast)\over k}.
```

Code location:

- `cpp/two_fluid_tt.cpp:977-989` evolves each `k` mode to recombination and
  stores the SW and Doppler source functions.

## ISW Term

The integrated Sachs-Wolfe contribution is computed as a separate line-of-sight
integral,

```math
\Theta_\ell^{\rm ISW}(k)
=
\int_{\eta_\ast}^{\eta_0}
2\Phi'(k,\eta)\,
j_\ell\!\left[k(\eta_0-\eta)\right]\,d\eta .
```

This term is not multiplied by the Silk or finite-width damping factors.

Code location:

- `cpp/two_fluid_tt.cpp:888-930` computes `Phi'(k, eta)` on the ISW time grid.
- `cpp/two_fluid_tt.cpp:932-974` performs the ISW line-of-sight integral.

## Projection Integral

For the primary source,

```math
\Theta_\ell^{\rm primary}(k)
=
\sqrt{D(k)}
\left[
S_{\rm SW}(k) j_\ell(k\chi_\ast)
+ S_{\rm Doppler}(k) j_\ell'(k\chi_\ast)
\right],
```

where

```math
\chi_\ast = \eta_0-\eta_\ast.
```

The damping factor is

```math
D(k)
=
\exp\left[-k^2\left(2x_s^2+c_s^2\sigma_\eta^2\right)\right].
```

Here `x_s^2 = k_D^{-2}` is computed from the photon diffusion integral

```math
k_D^{-2}(\eta_\ast)
=
\int_0^{\eta_\ast}
{d\eta\over 6\dot\kappa}
{R_b^2 + {16\over 15}(1+R_b)\over (1+R_b)^2},
```

where

```math
R_b={3\rho_b\over 4\rho_\gamma}
```

and

```math
\dot\kappa = a n_e \sigma_T .
```

The finite-width term uses the Gaussian width of the RECFAST visibility
function:

```math
\sigma_\eta^2
=
{\int d\eta\,g(\eta)(\eta-\bar\eta)^2
\over
\int d\eta\,g(\eta)}.
```

The full transfer function is

```math
\Theta_\ell(k)
=
\Theta_\ell^{\rm primary}(k)
+ \Theta_\ell^{\rm ISW}(k).
```

The TT power spectrum is computed as

```math
C_\ell^{TT}
=
4\pi A_s T_{\rm CMB}^2
\int {dk\over k}\,
\Theta_\ell^2(k).
```

The code outputs both `C_l` and

```math
{\ell(\ell+1)C_\ell\over 2\pi}
```

in units of microkelvin squared.

Code location:

- `cpp/two_fluid_tt.cpp:535-565` builds the source splines, `k` grid, ell
  grid, and normalization.
- `cpp/two_fluid_tt.cpp:566-578` computes the damping factor and samples the
  SW/Doppler source splines.
- `cpp/two_fluid_tt.cpp:580-603` prepares the spherical Bessel functions.
- `cpp/two_fluid_tt.cpp:612-642` evaluates the projection integral and the
  final `C_l`.
- `cpp/two_fluid_tt.cpp:1542-1558` interpolates CLASS-sampled output to every
  integer multipole in server mode.
- `cpp/two_fluid_tt.cpp:1560-1590` writes the spectrum to disk in command-line
  mode.

## Useful Runtime Options

The source components can be selected with

```text
--source-mode full
--source-mode sw
--source-mode doppler
--source-mode isw
```

The default ell sampling follows the CLASS-style sparse grid and then
interpolates to every integer ell in the output:

```text
--ell-grid class --interpolated-output
```

The default precision settings used for the current plots are

```text
--n-source 1001 --n-isw 320 --isw-k-eta0-max 900 --k-samples-per-period 40
```

The recombination redshift, finite-width damping scale, and Silk damping scale
are always computed from the RECFAST history.  There is no fixed or fitted
damping-scale option in this simplified package.

Normal command-line runs print a scale summary to stderr:

```text
z_star, eta_star, eta_0, chi_star, a_eq, k_eq, sigma_eta, k_width, k_D,
sound_speed_star
```

Server mode does not print this summary, so the Python wrapper protocol remains
simple.

Code location:

- `cpp/two_fluid_tt.cpp:26-63` defines the default parameters.
- `cpp/two_fluid_tt.cpp:507-545` computes and prints the scale summary.
- `cpp/two_fluid_tt.cpp:1465-1539` parses and checks command-line options.
