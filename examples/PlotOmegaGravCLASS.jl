using OmegaGrav
using MatterPower
using Plots, LaTeXStrings
using Dierckx
using DelimitedFiles, Printf

# %% Cosmological parameters
As, ns, kpivot = 2.097e-9, 0.9652, 0.05
Ωm, ΩΛ, Ωb = 0.315, 0.685, 0.049
Ωk = 1 - Ωm - ΩΛ
h0 = 0.674
ωm, fb = Ωm * h0^2, Ωb / Ωm

# %% Loop over redshifts
redshift = [0, 0.3, 0.5, 0.7, 1, 1.3, 1.5]
Ωgrav_pklin = zeros(7)
Ωgrav_pknl = zeros(7)
Ωgrav_halo = zeros(7)
for ired = 1:7
   # %% Read in linear and non-linear power spectra computed by CLASS
   filename = @sprintf("data/matterpower05_z%1d_pk.dat", ired)
   d = readdlm(filename, comments = true)
   pk = Spline1D(d[:, 1], d[:, 2])
   filename = @sprintf("data/matterpower05_z%1d_pk_nl.dat", ired)
   d = readdlm(filename, comments = true)
   pknl = Spline1D(d[:, 1], d[:, 2])

   # %% Compute Ωgrav
   z = redshift[ired]
   println("redshift = ", z)
   # Linear matter power spectrum
   Ωgrav_pklin[ired] = ograv_pk(x -> pk(x), z, Ωm)
   println("Ωgrav(linear pk) = ", Ωgrav_pklin[ired])
   # Non-linear matter power spectrum
   Ωgrav_pknl[ired] = ograv_pk(x -> pknl(x), z, Ωm)
   println("Ωgrav(non-linear pk) = ", Ωgrav_pknl[ired])
   # Halos
   Ωgrav_halo[ired] = ograv_halo(x -> pk(x), z, Ωm)
   println("Ωgrav(halos) = ", Ωgrav_halo[ired])
end

# %% Plot results and save to "omegagrav_class.pdf"
p = plot(
   redshift,
   -Ωgrav_pklin,
   xlab = "Redshift",
   ylab = L"-\Omega_{\rm grav}(z)~from~CLASS~P(k)",
   lab = "Linear P(k)",
   yaxis = :log10,
   legend = :bottomleft,
   ls = :dash,
   m = 2,
)
p = plot!(redshift, -Ωgrav_pknl, m = 2, lab = "Non-linear P(k)")
p = plot!(redshift, -Ωgrav_halo, m = 2, lab = "Halos")
savefig("omegagrav_class.pdf")
display(p)
