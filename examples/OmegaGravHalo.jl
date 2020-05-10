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

# %% Compute Ωgrav from halos, using M200m and Mvir
z = redshift
println("redshift = ", z)
Ωgrav_halo = ograv_halo(pk, z, Ωm)
println("Ωgrav(Δm=200) = ", Ωgrav_halo)
Ωgrav_halo = ograv_halo(pk, z, Ωm, virial = true)
println("Ωgrav(Δvir) = ", Ωgrav_halo)
