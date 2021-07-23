# OmegaGrav

This package contains functions to compute gravitational binding energy associated with large-scale clustering of matter (the so-called *large-scale structure*) in the Universe.

If you came here to reproduce the results in the paper ([Chiang, Makiya, Komatsu & Ménard, arXiv:2007.01679](https://arxiv.org/abs/2007.01679)), you can use [examples/PlotFigure1.jl](https://github.com/komatsu5147/OmegaGrav.jl/tree/master/examples/PlotFigure1.jl), [examples/PlotFigure2.jl](https://github.com/komatsu5147/OmegaGrav.jl/tree/master/examples/PlotFigure2.jl),
[examples/PlotFigure3.jl](https://github.com/komatsu5147/OmegaGrav.jl/tree/master/examples/PlotFigure3.jl), [examples/PlotFigure4.jl](https://github.com/komatsu5147/OmegaGrav.jl/tree/master/examples/PlotFigure4.jl), [examples/GenerateTable1.jl](https://github.com/komatsu5147/OmegaGrav.jl/tree/master/examples/GenerateTable1.jl), and [examples/GenerateTable2.jl](https://github.com/komatsu5147/OmegaGrav.jl/tree/master/examples/GenerateTable2.jl) to reproduce Figure 1-4 and Table 1,2 of the paper, respectively. Running these codes requires the python wrapper for CLASS [classy](https://github.com/lesgourg/class_public/wiki/Python-wrapper). The script [examples/compute_pk_class.jl](https://github.com/komatsu5147/OmegaGrav.jl/tree/master/examples/compute_pk_class.jl) calls `classy` via `PyCall`.

To calculate the Compton-y-weighted halo clustering bias, `b_y`, given in Equation (B1) of [Chiang, Makiya, Ménard & Komatsu, arXiv:2006.14650](https://arxiv.org/abs/2006.14650), you can use [examples/GenerateBy.jl](https://github.com/komatsu5147/OmegaGrav.jl/tree/master/examples/GenerateBy.jl).

## Installation

From the Julia REPL, run
```Julia
using Pkg
Pkg.add(url="https://github.com/komatsu5147/OmegaGrav.jl")
```
## Contents

The package contains
- ``ograv_pk(pk, z, Ωm; kmin=5e-4, kmax=3e1)``: returns the comoving density parameter, `Ωgrav = Ωm * W / 2`, of gravitational binding energy computed from a matter power spectrum `pk`. Here, `W` is the mean gravitational potential energy per unit mass (see "The physics behind" at the bottom of this page). It is based on Equation (60) of [Fukugita & Peebles, ApJ, 616, 643 (2004)](https://iopscience.iop.org/article/10.1086/425155), extended to arbitrary redshift `z` by Equation (14) of [Chiang, Makiya, Komatsu & Ménard, arXiv:2007.01679](https://arxiv.org/abs/2007.01679).
- ``ograv_halo(pk, z, Ωm[, Ωcb=Ωm]; Mmin=1e11, Mmax=5e15, virial=false, t10MF=false)``: returns `Ωgrav` from gravitationally collapsed structures (*halos*). It uses Equation (17) of [Chiang, Makiya, Komatsu & Ménard, arXiv:2007.01679](https://arxiv.org/abs/2007.01679).

## Arguments
- `pk::Any`(k): a function which returns a linear matter power spectrum with the argument k being the comoving wavenumber. This can be an interpolation function constructed from tabulated data.
   - `pk` is in units of (Mpc/h)^3 and the argument k is in units of h/Mpc.
- `z::Real`: redshift.
- `Ωm::Real`: present-day total matter density parameter.

## Optional arguments
- `Ωcb::Real`: present-day baryon + cold dark matter density parameter. The default value is equal to `Ωm` given in the argument. This parameter is relevant when the neutrino contribution to the mass density needs to be excluded.

## Optional keyword arguments
- `kmin::Real`: minimum wavenumber for integration, ``∫_{kmin}^{kmax} dk P(k)``. The default value: `5e-4`.
- `kmax::Real`: maximum wavenumber for integration, ``∫_{kmin}^{kmax} dk P(k)``. The default value: `3e1`.
- `Mmin::Real`: minimum mass for integration, ``∫_{Mmin}^{Mmax} dM dn/dM Ag GM^2/R``. The default value: `1e11`.
- `Mmax::Real`: maximum mass for integration, ``∫_{Mmin}^{Mmax} dM dn/dM Ag GM^2/R``. The default value: `5e15`.
- `virial::Bool`: if `true`, use the virial overdensity `Δvir`. If `false` (the default), use `Δm=200`.
- `t10MF::Bool`: if `true`, use `tinker10MF` for the halo multiplicity function. If `false` (the default), use `tinker08MF`. See [HaloMF.jl](https://github.com/komatsu5147/HaloMF.jl) for definition of these functions.

## Example Juia codes
### 1. Ωgrav from linear and non-linear matter power spectra

The following example code (avaiable in [examples/OmegaGravPk.jl](https://github.com/komatsu5147/OmegaGrav.jl/blob/master/examples/OmegaGravPk.jl)) computes Ωgrav from linear and non-linear power spectra.

If you would like to generate a nice figure showing Ωgrav as a function of redshift, take a look at [examples/PlotOmegaGrav.jl](https://github.com/komatsu5147/OmegaGrav.jl/blob/master/examples/PlotOmegaGrav.jl).

These codes use an approximate "no-wiggle" linear matter transfer function (see [MatterPower.jl](https://github.com/komatsu5147/MatterPower.jl) for descripton). This is sufficiently accurate for most purposes, but we also provide example codes to calculate and plot Ωgrav from tabulated data of the power spectrum generated by [CLASS](https://github.com/lesgourg/class_public). See [examples/OmegaGravPkCLASS.jl](https://github.com/komatsu5147/OmegaGrav.jl/blob/master/examples/OmegaGravPkCLASS.jl) and [examples/PlotOmegaGravCLASS.jl](https://github.com/komatsu5147/OmegaGrav.jl/blob/master/examples/PlotOmegaGravCLASS.jl)
```Julia
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

```Julia
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
## The physics behind

What is *binding energy*? Consider a hydrogen atom: it has a proton and an electron. This is a bound system, and it takes some energy to remove the electron from the atom. This energy is called binding energy. For a hydrogen atom, the binding energy is -13.6 electron volts. Note the negative sign; you need to inject energy to remove the electron. When the injected energy and binding energy sum to zero, the electron is able to escape the atom. Any extra energy above that will go to kinetic energy of the escaping electron.

The same applies to a gravitationally bound system. Consider a sphere with a radius R and a uniform density ρ. It is known (see [wikipedia](https://en.wikipedia.org/wiki/Gravitational_binding_energy)) that gravitational binding energy of a uniform density sphere is -3GM<sup>2</sup>/(5R). This is the energy required to tear a gravitating uniform density sphere apart.

We can define the total potential energy associated with an arbitrary distribution of matter in expanding space. Let us write the total potential energy per unit mass as W and the total mass of a system as M. Equation (24.2) of Jim Peeble's 1980 book "[Large-scale Structure of the Universe](https://press.princeton.edu/books/paperback/9780691209838/the-large-scale-structure-of-the-universe)" gives

MW = (1/2)a<sup>3</sup>∫ d<sup>3</sup>x ϕ(x) [ρ(x) - ρ<sub>b</sub>] -- (1)

where x denotes the comoving coordinates, a is the scale factor, a<sup>3</sup>d<sup>3</sup>x is a physical volume element, ϕ(x) and ρ(x) are distributions of Newton's gravitational potential and matter density, respectively, and ρ<sub>b</sub> is the background matter density. Thus, the potential is sourced by a density excess above the background, δρ ≡ ρ(x) - ρ<sub>b</sub>. For a uniform density sphere above the background with a physical radius R,
a mass excess δM = 4πG δρ R<sup>3</sup>/3 and ϕ(r≦R) = -2πG δρ (R<sup>2</sup>-r<sup>2</sup>/3), we find MW = -3G(δM)<sup>2</sup>/(5R). This reproduces the factor of 3/5 in gravitaional binding energy of a uniform density sphere.

The potential ϕ is related to the density excess via Poisson equation:

∇<sup>2</sup>ϕ = 4πG a<sup>2</sup> δρ -- (2)

where a spatial differential operator ∇ is defined with respect to the comoving coordinates x (hence the factor a<sup>2</sup> in the right hand side).
Now, Fourier transforming Equations (1,2), taking ensemble average of Eq.(1) and dividing both sides
of Eq.(1) by the total mass of the system M = ρ<sub>b</sub> a<sup>3</sup>∫d<sup>3</sup>x, we obtain

W = -(1/π)G a<sup>2</sup>ρ<sub>b</sub> ∫dk P(k) -- (3)

where P(k) is the power spectrum of excess matter density constast, δρ/ρ<sub>b</sub>. Eq.(3) is Fourier transform of Equation (24.9) of  [Peeble's book](https://press.princeton.edu/books/paperback/9780691209838/the-large-scale-structure-of-the-universe). Thus, once we have P(k) (from, e.g., [MatterPower.jl](https://github.com/komatsu5147/MatterPower.jl)), we can integrate Eq.(3) to find the total gravitational potential energy of the large-scale structure in the Universe. Isn't it neat?

Eq.(3) receives contributions from excess matter clustering at all scales. However, not all structures in the Universe have collapsed gravitationally, like galaxies and clusters of galaxies. We name gravitationally collapsed structures in the Universe collectively as "*halos*". How can we isolate contributions from halos? The matter power spectrum can be divided into contributions from matter clustering inside halos ("1-halo term", P<sub>1h</sub>) and clustering between different halos ("2-halo term", P<sub>2h</sub>). The 1-halo term is given by

P<sub>1h</sub>(k) = (1/ρ<sub>b0</sub><sup>2</sup>)∫ dM dn/dM M<sup>2</sup> |u(k,R)|<sup>2</sup>  -- (4)

where k is the comoving wavenumber, ρ<sub>b0</sub> = a<sup>3</sup>ρ<sub>b</sub> is the present-day background matter density, M is the mass of a halo, dn/dM is the number density of halos per unit mass interval (see [HaloMF.jl](https://github.com/komatsu5147/HaloMF.jl)), and u(k,M) is Fourier transform of density profile of halos normalized as u(k,M) → 1 for k → 0.  

Using Eq.(4) in Eq.(3), we find

W<sub>1h</sub> = -(1/ρ<sub>b0</sub>)∫ dM dn/dM  (GM<sup>2</sup>/R) A<sub>g</sub>(M) -- (5)

where R is a physical size of a halo and A<sub>g</sub>(M) ≡ (1/π)∫ dk |u(k,M)|<sup>2</sup>. A<sub>g</sub> is equal to 3/5 for a uniform density sphere  and is independent of M, but it depends weakly on M in general. Equation (5) was derived by [Chiang, Makiya, Komatsu & Ménard, arXiv:2007.01679](https://arxiv.org/abs/2007.01679).

In this Julia package `OmegaGrav.jl`, we provide functions `ograv_pk(pk, z, Ωm)` and `ograv_halo(pk, z, Ωm)` to evaluate Equations (3) and (5), respectively, as a function of redshift `z`. More precisely, the functions return the comoving density parameter `Ωgrav` defined by ``Ωgrav = Ωm W/2`` (following Section 2.4 of [Fukugita & Peebles, ApJ, 616, 643 (2004)](https://iopscience.iop.org/article/10.1086/425155)), where `Ωm` is the present-day matter density parameter. For example, the calculation using these Julia functions shows that `ograv_halo` is about 1/3 and 1/10 of `ograv_pk` with the non-linear P(k) at z=0 and 1.5, respectively. Good to know!

## Acknowledgment

Development of the functions provided in this package was supported in part by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy - EXC-2094 - 390783311 and JSPS KAKENHI Grant Number JP15H05896.
