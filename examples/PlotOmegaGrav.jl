using OmegaGrav
using MatterPower
using Plots, LaTeXStrings

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
   # %% Tabulate linear growth factor as a function of scale factor, a
   sol = setup_growth(Ωm, ΩΛ)
   a = 1 / (1 + redshift[ired])
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
   z = redshift[ired]
   E2 = Ωm * (1 + z)^3 + Ωk * (1 + z)^2 + ΩΛ
   Ωmz = Ωm * (1 + z)^3 / E2
   # Define a function to return a halofit non-linear power spectrum
   pknl(kovh) = halofit(pk, p, Ωmz, kovh) # Mpc^3/h^3

   # %% Compute Ωgrav
   println("redshift = ", z)
   # Linear matter power spectrum
   Ωgrav_pklin[ired] = ograv_pk(pk, z, Ωm)
   println("Ωgrav(linear pk) = ", Ωgrav_pklin[ired])
   # Non-linear matter power spectrum
   Ωgrav_pknl[ired] = ograv_pk(pknl, z, Ωm)
   println("Ωgrav(non-linear pk) = ", Ωgrav_pknl[ired])
   # Halos
   Ωgrav_halo[ired] = ograv_halo(pk, z, Ωm)
   println("Ωgrav(halos) = ", Ωgrav_halo[ired])
end

# %% Plot results and save to "omegagrav.pdf"
p = plot(
   redshift,
   -Ωgrav_pklin,
   xlab = "Redshift",
   ylab = L"-\Omega_{\rm grav}(z)",
   lab = "Linear P(k)",
   yaxis = :log10,
   legend = :bottomleft,
   ls = :dash,
   m = 2,
)
p = plot!(redshift, -Ωgrav_pknl, m = 2, lab = "Non-linear P(k)")
p = plot!(redshift, -Ωgrav_halo, m = 2, lab = "Halos")
savefig("omegagrav.pdf")
display(p)
