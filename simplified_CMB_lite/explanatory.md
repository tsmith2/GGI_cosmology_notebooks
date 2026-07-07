# simplified_CMB_lite Equations

This lite version keeps the same RECFAST-based recombination and damping
calculation as `simplified_CMB`, but removes the ISW line-of-sight integral.
The temperature projection uses only the instantaneous last-scattering
Sachs-Wolfe and Doppler sources.

## Background

The code evolves the scale factor in conformal time:

```math
a'(\eta) =
a^2 H_0
\left[
\Omega_m a^{-3} + \Omega_\gamma a^{-4} + \Omega_\Lambda
\right]^{1/2}.
```

The model is spatially flat:

```math
\Omega_\Lambda = 1 - \Omega_m - \Omega_\gamma .
```

Implementation: see `background_rhs`, `compute_background_solution`, and
`Background`.

## Recombination And Damping

The compact RECFAST-style module computes `x_e(z)`.  From that history the code
finds the visibility peak and derives:

```math
z_\star,\qquad \eta_\star,\qquad
\sigma_\eta,\qquad k_D .
```

The primary source is damped by

```math
\exp\left[-k^2\left(2 k_D^{-2}+c_s^2\sigma_\eta^2\right)\right].
```

Implementation: see `compute_recfast_history` and
`compute_recombination_scales_from_recfast`.

## Perturbations

The dynamical variables are

```math
\left(\delta_\gamma,\delta_c,\theta_{\gamma b},\theta_c,\Phi\right).
```

During tight coupling,

```math
\delta_b = {3\over4}\delta_\gamma .
```

The equations are

```math
\delta_\gamma' = -{4\over3}\theta_{\gamma b}+4\Phi',
```

```math
\delta_c' = -\theta_c+3\Phi',
```

```math
\theta_{\gamma b}' =
-{\mathcal H\over 1+R}\theta_{\gamma b}
+ k^2\Phi
+ {R\over4(1+R)}k^2\delta_\gamma,
```

```math
\theta_c' = -\mathcal H\theta_c+k^2\Phi,
```

with

```math
R={4\rho_\gamma\over3\rho_b},\qquad
\mathcal H={a'\over a}.
```

The Newtonian potential is obtained from

```math
k^2\Phi + 3\mathcal H(\Phi' + \mathcal H\Phi)
= -{3\over2}\mathcal H^2\delta_{\rm tot},
```

where

```math
\delta_{\rm tot}
={\rho_\gamma\delta_\gamma+\rho_b\delta_b+\rho_c\delta_c\over
\rho_\gamma+\rho_b+\rho_c+\rho_\Lambda}.
```

Implementation: see `rhs`, `initial_conditions`, and `integrate_source`.

## Initial Conditions

The default superhorizon radiation-era initial conditions are

```math
\Phi_i=-{2\over3},\qquad
\delta_{\gamma,i}=-2\Phi_i,\qquad
\delta_{c,i}=\delta_{b,i}={3\over4}\delta_{\gamma,i},
```

```math
\theta_{\gamma b,i}=\theta_{c,i}
={1\over2} k^2\eta_i\Phi_i .
```

Here `Phi_i = -2/3` is chosen so that, on superhorizon scales during radiation
domination, the potential normalization corresponds to the primordial curvature
normalization used by the teaching code.

## Temperature Projection

For each wavenumber, the source at recombination is

```math
S_{\rm SW}(k)=\Phi(\eta_\star,k)+{1\over4}\delta_\gamma(\eta_\star,k),
```

```math
S_{\rm Dop}(k)=-{\theta_{\gamma b}(\eta_\star,k)\over k}.
```

The transfer function is

```math
\Theta_\ell(k)
= D(k)\left[S_{\rm SW}(k)j_\ell(k\chi_\star)
+S_{\rm Dop}(k)j_\ell'(k\chi_\star)\right],
```

with

```math
\chi_\star=\eta_0-\eta_\star .
```

The TT spectrum is

```math
C_\ell^{TT}
= 4\pi A_s T_{\rm CMB}^2
\int {dk\over k}\,\Theta_\ell^2(k).
```

Implementation: see `solve_sources` and `compute_spectrum`.

## Speed Choices

Mode evolution is parallelized over independent wavenumbers in `solve_sources`.
Use `--n-threads 0` for automatic threading, or set a positive integer to force
a specific thread count.

The default fast path uses a fixed-argument spherical Bessel cache through
`--bessel-x-cache`.  The first run creates the cache; later runs load it.
