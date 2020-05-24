using OmegaGrav
using MatterPower
using Dierckx
using DelimitedFiles, Printf

# %% Cosmological parameters
Ωm = 0.315

# %% Specify a redshift
redshift = [0, 0.3, 0.5, 0.7, 1, 1.3, 1.5]
ired = 1
z = redshift[ired]

# %% Read in linear and non-linear power spectra computed by CLASS
filename = @sprintf("data/matterpower05_z%1d_pk.dat", ired)
d = readdlm(filename, comments = true)
pk = Spline1D(d[:, 1], d[:, 2])
filename = @sprintf("data/matterpower05_z%1d_pk_nl.dat", ired)
d = readdlm(filename, comments = true)
pknl = Spline1D(d[:, 1], d[:, 2])

# %% Compute Ωgrav
println("redshift = ", z)
# Linear matter power spectrum
Ωgrav_pklin = ograv_pk(x -> pk(x), z, Ωm)
println("Ωgrav(linear pk) = ", Ωgrav_pklin)
# Non-linear matter power spectrum
Ωgrav_pknl = ograv_pk(x -> pknl(x), z, Ωm)
println("Ωgrav(non-linear pk) = ", Ωgrav_pknl)
