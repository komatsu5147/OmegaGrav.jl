# OmegaGrav

This package contains functions to compute gravitational binding energy associated with large-scale clustering of matter (the so-called *large-scale structure*) in the Universe.

The package contains
- ``ograv_pk(pk, z, Ωm; kmin=5e-4, kmax=1e2)``: returns the comoving density parameter, `Ωgrav`, of gravitational binding energy computed from a matter power spectrum `pk`. It is based on Equation (60) of [Fukugita & Peebles, ApJ, 616, 643 (2004)](https://iopscience.iop.org/article/10.1086/425155), extended to arbitrary redshift `z` by Chiang, Makiya, Komatsu & Ménard (in prep)
- ``ograv_halo(pk, z, Ωm; Mmin=5e8, Mmax=5e15, virial=false, t10MF=false)``: returns `Ωgrav` from gravitationally collapsed structures (*halos*). It uses Equation (TBD) of Chiang, Makiya, Komatsu & Ménard (in prep)

## Arguments
- `pk::Any`(k): a function which returns a linear matter power spectrum with the argument k being the comoving wavenumber. This can be an interpolation function constructed from tabulated data.
- `z::Real`: redshift.
- `Ωm::Real`: present-day total matter density parameter.

## Optional arguments
- `kmin::Real`: minimum wavenumber for integration, ``∫_{kmin}^{kmax} dk P(k)``. The default value: `5e-4`.
- `kmax::Real`: maximum wavenumber for integration, ``∫_{kmin}^{kmax} dk P(k)``. The default value: `1e2`.
- `Mmin::Real`: minimum mass for integration, ``∫_{Mmin}^{Mmax} dM dn/dM Ag GM^2/R``. The default value: `5e8`.
- `Mmax::Real`: maximum mass for integration, ``∫_{Mmin}^{Mmax} dM dn/dM Ag GM^2/R``. The default value: `5e15`.
- `virial::Bool=false`: if `true`, use the virial overdensity `Δvir`. If `false` (the default), use `Δm=200`.
- `t10MF::Bool=false`: if `true`, use `tinker10MF` for the halo multiplicity function. If `false` (the default), use `tinker08MF`. See [HaloMF.jl](https://github.com/komatsu5147/HaloMF.jl) for definition of these functions.

## Example Juia codes
### 1. Ωgrav from linear and non-linear matter power spectra

The following example code (avaiable in [examples/OmegaGravPk.jl](https://github.com/komatsu5147/OmegaGrav.jl/blob/master/examples/OmegaGravPk.jl)) computes Ωgrav from linear and non-linear power spectra.

If you would like to generate a nice figure showing Ωgrav as a function of redshift, take a look at [examples/PlotOmegaGrav.jl](https://github.com/komatsu5147/OmegaGrav.jl/blob/master/examples/PlotOmegaGrav.jl).

These codes use an approximate "no-wiggle" linear matter transfer function (see [MatterPower.jl](https://github.com/komatsu5147/MatterPower.jl) for descripton). This is sufficiently accurate for most purposes, but we also provide codes to read in power spectra calculated by [CLASS](https://github.com/lesgourg/class_public). See [examples/OmegaGravPkCLASS.jl](https://github.com/komatsu5147/OmegaGrav.jl/blob/master/examples/OmegaGravPkCLASS.jl) and [examples/PlotOmegaGravCLASS.jl](https://github.com/komatsu5147/OmegaGrav.jl/blob/master/examples/PlotOmegaGravCLASS.jl)
```
using OmegaGrav
using MatterPower

# %% Specify a redshift
redshift = 0

# %% Cosmological parameters
As, ns, kpivot = 2.097e-9, 0.9652, 0.05
Ωm, ΩΛ, Ωb = 0.315, 0.685, 0.049
Ωk = 1 - Ωm - ΩΛ
h0 = 0.674
ωm, fb = Ωm * h0^2, Ωb / Ωm

# %% Tabulate linear growth factor as a function of scale factor, a
sol = setup_growth(Ωm, ΩΛ)
a = 1 / (1 + redshift)
D1 = sol(a)[1]

# %% Define a function to return a linear matter power spectrum (in units of Mpc^3/h^3)
# as a function of the comoving wavenumber, kovh, in units of h/Mpc.
# Here is an example using Einstein & Hu's analytical transfer function
pk(kovh) =
   D1^2 *
   As *
   (kovh * h0 / kpivot)^(ns - 1) *
   (2 * kovh^2 * 2998^2 / 5 / Ωm)^2 *
   t_nowiggle(kovh * h0, ωm, fb)^2 *
   2 *
   π^2 / kovh^3

# %% Compute the non-linear power spectrum using halofit
# Get three parameters [kσ, neff, C] needed for the halofit transform
p = kσ, neff, C = setup_halofit(pk)
# Calculate Ωm(z) = Ωm(z=0)(1+z)^3 / E^2(z), where E^2(z) = H^2(z) / H0^2
z = redshift
E2 = Ωm * (1 + z)^3 + Ωk * (1 + z)^2 + ΩΛ
Ωmz = Ωm * (1 + z)^3 / E2
# Define a function to return a halofit non-linear power spectrum
pknl(kovh) = halofit(pk, p, Ωmz, kovh) # Mpc^3/h^3

# %% Compute Ωgrav
println("redshift = ", z)
# Linear matter power spectrum
Ωgrav_pklin = ograv_pk(pk, z, Ωm)
println("Ωgrav(linear pk) = ", Ωgrav_pklin)
# Non-linear matter power spectrum
Ωgrav_pknl = ograv_pk(pknl, z, Ωm)
println("Ωgrav(non-linear pk) = ", Ωgrav_pknl)
```
### 2. Ωgrav from halos

The following example code (avaiable in [examples/OmegaGravHalo.jl](https://github.com/komatsu5147/OmegaGrav.jl/blob/master/examples/OmegaGravHalo.jl)) computes Ωgrav from halos.

To test robustness of the results against various options, the code computes Ωgrav with
- `tinker08MF` and overdensity `Δm = 200`
- `tinker08MF` and overdensity `Δvir`
- `tinker10MF` and overdensity `Δm = 200`
- `tinker10MF` and overdensity `Δvir`

```
using OmegaGrav
using MatterPower

# %% Specify a redshift
redshift = 0

# %% Cosmological parameters
As, ns, kpivot = 2.097e-9, 0.9652, 0.05
Ωm, ΩΛ, Ωb = 0.315, 0.685, 0.049
Ωk = 1 - Ωm - ΩΛ
h0 = 0.674
ωm, fb = Ωm * h0^2, Ωb / Ωm

# %% Tabulate linear growth factor as a function of scale factor, a
sol = setup_growth(Ωm, ΩΛ)
a = 1 / (1 + redshift)
D1 = sol(a)[1]

# %% Define a function to return a linear matter power spectrum (in units of Mpc^3/h^3)
# as a function of the comoving wavenumber, kovh, in units of h/Mpc.
# Here is an example using Einstein & Hu's analytical transfer function
pk(kovh) =
   D1^2 *
   As *
   (kovh * h0 / kpivot)^(ns - 1) *
   (2 * kovh^2 * 2998^2 / 5 / Ωm)^2 *
   t_nowiggle(kovh * h0, ωm, fb)^2 *
   2 *
   π^2 / kovh^3

# %% Compute Ωgrav from halos, using M200m and Mvir, and tinker08MF and tinker10MF
z = redshift
println("redshift = ", z)
Ωgrav_halo = ograv_halo(pk, z, Ωm)
println("Ωgrav(Δm=200, T08) = ", Ωgrav_halo)
Ωgrav_halo = ograv_halo(pk, z, Ωm, virial = true)
println("Ωgrav(Δvir, T08) = ", Ωgrav_halo)
Ωgrav_halo = ograv_halo(pk, z, Ωm, t10MF = true)
println("Ωgrav(Δm=200, T10) = ", Ωgrav_halo)
Ωgrav_halo = ograv_halo(pk, z, Ωm, virial = true, t10MF = true)
println("Ωgrav(Δvir, T10) = ", Ωgrav_halo)
```
