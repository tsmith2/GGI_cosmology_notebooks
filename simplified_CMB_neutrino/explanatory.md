# Simplified CMB With Perfect-Fluid Neutrinos

This package evolves a simplified flat LCDM photon-baryon-CDM-neutrino model
in Newtonian gauge.  The independent time variable is conformal time `eta`, in
Mpc, and the wavenumber is the physical `k`, in Mpc^-1.

The evolved perturbation state is

```text
u = (delta_gamma, delta_nu, delta_c, theta_gammab, theta_nu, theta_c, Phi).
```

The baryon density perturbation is not independently evolved.  It is set by
the tight-coupling relation

```math
\delta_b = {3\over 4}\delta_\gamma .
```

Code location:

- `cpp/two_fluid_tt.cpp:26-60` defines the default parameters, including
  `N_eff = 3`.
- `cpp/two_fluid_tt.cpp:62-77` computes `Omega_nu`, total radiation, and
  flatness.
- `cpp/two_fluid_tt.cpp:233-234` defines the seven-component perturbation
  state.

## Background

The massless-neutrino density is

```math
\Omega_\nu =
{7\over 8}\left({4\over 11}\right)^{4/3}
N_{\rm eff}\Omega_\gamma .
```

The radiation density used in the background is

```math
\Omega_r=\Omega_\gamma+\Omega_\nu .
```

The model is strictly flat:

```math
\Omega_\Lambda = 1-\Omega_c-\Omega_b-\Omega_r .
```

The conformal Hubble parameter is

```math
{\cal H}^2
=
a^2 H_0^2
\left(
{\Omega_r\over a^4}
+{\Omega_b\over a^3}
+{\Omega_c\over a^3}
+\Omega_\Lambda
\right).
```

The conformal time is integrated from

```math
{d\eta\over d\ln a}={1\over aH(a)}.
```

Code location:

- `cpp/two_fluid_tt.cpp:130-142` uses total radiation in the RECFAST
  expansion rate.
- `cpp/two_fluid_tt.cpp:368-429` integrates the background ODE.
- `cpp/two_fluid_tt.cpp:476-488` computes the RECFAST history and derives
  `z_star`, `sigma_eta`, and the Silk scale.

## Perturbations

The total density contrast entering the Einstein constraint is

```math
\delta_{\rm tot}
=
{\rho_\gamma\delta_\gamma
+\rho_\nu\delta_\nu
+\rho_b\delta_b
+\rho_c\delta_c
\over
\rho_\gamma+\rho_\nu+\rho_b+\rho_c+\rho_\Lambda}.
```

The potential is evolved by solving

```math
k^2\Phi+3{\cal H}\left(\Phi'+{\cal H}\Phi\right)
=
-{3\over 2}{\cal H}^2\delta_{\rm tot}
```

for `Phi'`.

The photon-baryon inertia ratio is

```math
R={4\rho_\gamma\over 3\rho_b}.
```

The evolved perturbation equations are

```math
\delta_\gamma'=-{4\over 3}\theta_{\gamma b}+4\Phi',
```

```math
\delta_\nu'=-{4\over 3}\theta_\nu+4\Phi',
```

```math
\delta_c'=-\theta_c+3\Phi',
```

```math
\theta_{\gamma b}'
=
-{{\cal H}\over 1+R}\theta_{\gamma b}
+k^2\Phi
+{R\over 4(1+R)}k^2\delta_\gamma,
```

```math
\theta_\nu'=k^2\Phi+{1\over 4}k^2\delta_\nu,
```

```math
\theta_c'=-{\cal H}\theta_c+k^2\Phi.
```

Code location:

- `cpp/two_fluid_tt.cpp:752-779` evaluates the perturbation RHS.

## Initial Conditions

The initial time is

```math
\eta_i=\min\left({10^{-3}\over k},10^{-3}\eta_{\rm eq}\right).
```

The initial conditions are

```math
\Phi_i=-{2\over 3},
```

```math
\delta_{\gamma,i}=\delta_{\nu,i}=-2\Phi_i={4\over 3},
```

```math
\delta_{b,i}=\delta_{c,i}
={3\over 4}\delta_{\gamma,i}=1,
```

```math
\theta_{\gamma b,i}=\theta_{\nu,i}=\theta_{c,i}
={1\over 2}k^2\eta_i\Phi_i
=-{1\over 3}k^2\eta_i.
```

Code location:

- `cpp/two_fluid_tt.cpp:792-799` sets the initial time and state.

## Projection

The instantaneous recombination Sachs-Wolfe and Doppler sources are

```math
S_{\rm SW}(k)=\Phi(k,\eta_\ast)+{1\over 4}\delta_\gamma(k,\eta_\ast),
```

```math
S_{\rm Doppler}(k)=-{\theta_{\gamma b}(k,\eta_\ast)\over k}.
```

The ISW source is

```math
\Theta_\ell^{\rm ISW}(k)
=
\int_{\eta_\ast}^{\eta_0}
2\Phi'(k,\eta)
j_\ell[k(\eta_0-\eta)]\,d\eta .
```

The full transfer function is

```math
\Theta_\ell(k)
=
\sqrt{D(k)}
\left[
S_{\rm SW}j_\ell(k\chi_\ast)
+S_{\rm Doppler}j_\ell'(k\chi_\ast)
\right]
+\Theta_\ell^{\rm ISW}(k),
```

with

```math
D(k)=\exp[-k^2(2x_s^2+c_s^2\sigma_\eta^2)].
```

The TT spectrum is

```math
C_\ell^{TT}
=
4\pi A_s T_{\rm CMB}^2
\int {dk\over k}\Theta_\ell^2(k).
```

Code location:

- `cpp/two_fluid_tt.cpp:535-642` computes the projection integral.
- `cpp/two_fluid_tt.cpp:900-984` computes the direct ISW line-of-sight term.
- `cpp/two_fluid_tt.cpp:987-999` computes the SW and Doppler source functions.

## Runtime Notes

The default neutrino parameter is

```text
--N-eff 3
```

The option can also be written as `--n-eff`.

Normal command-line runs print `N_eff`, `Omega_nu`, `z_star`, `k_eq`, `k_D`,
`k_width`, and related background scales.  Server mode stays quiet so the
Python wrapper protocol remains simple.
