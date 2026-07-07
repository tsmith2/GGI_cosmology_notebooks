# simplified_CMB_neutrino

This is a simplified RECFAST-enabled CMB temperature solver using
photon-CDM-neutrino theta-variable perturbation equations.  Neutrinos are
treated as a massless perfect fluid with default `N_eff = 3`.  Baryon density
perturbations are kept in the gravitational source, but are not independently
evolved:

```math
\delta_b = {3\over 4}\delta_\gamma .
```

This is approximately true during tight coupling.  The included compact
RECFAST-style module computes the recombination history.  The CMB temperature
projection still assumes instantaneous recombination, but the code now uses
RECFAST to compute `z_star`, the finite-width damping scale, and the Silk
damping scale.

The neutrino density is

```math
\Omega_\nu =
{7\over 8}\left({4\over 11}\right)^{4/3}N_{\rm eff}\Omega_\gamma .
```

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
  --N-eff 3 --bessel-class-memory --output smoke_tt.dat
```

See `explanatory.md` for the physical equations and corresponding C++ line
numbers.
