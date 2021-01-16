using OmegaGrav
using CSV, DataFrames
using Dierckx
using Plots, Plots.PlotMeasures, LaTeXStrings
# %% Compute the matter power spectrum using CLASS
include("compute_pk_class.jl")
# %% Compute `Ωgrav = Ωm * W / 2` and Ωtherm at redshifts of the data points
d = CSV.read("data/d16_Omega_th_data.csv", DataFrame)
nred = size(d)[1]
redshift = zeros(nred + 1)
Ωghalo = zeros(nred + 1)
Ωgpk = zeros(nred + 1)
Ωglin = zeros(nred + 1)
Ωth = zeros(nred + 1)
Ωthu = zeros(nred + 1)
Ωthl = zeros(nred + 1)
for ired = 1:nred+1
   if ired == 1
      z = 0 # Add z=0 which is not included in the data file
   else
      z = d.z[ired-1]
   end
   redshift[ired] = z
   # Define functions to return linear and non-linear power spectra
   # Note: The CLASS code takes wavenumbers in units of 1/Mpc (no h) and
   # return power spectra in units of Mpc^3 (no 1/h^3).
   pkcb_class(kovh) = cosmo.pk_cb_lin(kovh * h0, z) * h0^3
   pknl_class(kovh) = cosmo.pk(kovh * h0, z) * h0^3
   pklin_class(kovh) = cosmo.pk_lin(kovh * h0, z) * h0^3
   # Spline interpolate in log(k)
   lnk = log(1e-4):0.05:log(100)
   pkcb = Spline1D(lnk, pkcb_class.(exp.(lnk)))
   pknl = Spline1D(lnk, pknl_class.(exp.(lnk)))
   pklin = Spline1D(lnk, pklin_class.(exp.(lnk)))
   # %% Compute Ωgrav from non-linear total matter P(k)
   Ωgpk[ired] = ograv_pk(x -> pknl(log(x)), z, Ωm)
   # %% Compute Ωgrav from linear total matter P(k)
   Ωglin[ired] = ograv_pk(x -> pklin(log(x)), z, Ωm)
   # %% Compute Ωgrav from Halos, excluding the neutrino contribution
   Ωghalo[ired] = ograv_halo(x -> pkcb(log(x)), z, Ωm, Ωcb)
   # %% Compute Ωtherm from Halos, excluding the neutrino contribution
   Ωth[ired] = otherm_upp(x -> pkcb(log(x)), z, Ωm, h0, Ωcb)
   # %% Compute Ωtherm for upper and lower 68% confidence level in B
   Ωthl[ired] = otherm_upp(x -> pkcb(log(x)), z, Ωm, h0, Ωcb, massbias = 1.315)
   Ωthu[ired] = otherm_upp(x -> pkcb(log(x)), z, Ωm, h0, Ωcb, massbias = 1.221)
end

#%% Plot results and save to figure1.pdf
fb = params["omega_b"] / (params["omega_cdm"] + params["omega_b"])
Ωb = params["omega_b"] / h0^2
ΔΩg = Ωgpk - Ωglin
y = 0.5 * (Ωghalo + ΔΩg) * (-2fb / 3)
Δyu = ΔΩg * (-2fb / 3) - y
Δyl = y - Ωghalo * (-2fb / 3)
p = plot(
   redshift,
   -2 * Ωgpk,
   c = :black,
   lab = L"-\Omega_W",
   lw = 2,
   xlims = [0, 1.5],
   ylims = [1e-9, 2e-6],
   yaxis = :log,
   xlab = "Redshift, z",
   ylab = "Comoving Energy Density Parameters",
   legend = :topright,
   legendfontsize = 12,
   labelfontsize = 14,
   right_margin = 45px,
)
p = plot!(
   redshift,
   -2 * Ωghalo,
   c = :black,
   ls = :dash,
   lab = L"-\Omega_W^{halo}",
   lw = 2,
)
p = plot!(redshift, y, ribbon = (Δyl, Δyu), lab = L"\Omega_W^{ref}", c = :green)
p = scatter!(
   (d.z, d.Omega_th),
   yerror = (d.Omega_th .- d.Omega_th_low, d.Omega_th_up .- d.Omega_th),
   ms = 5,
   c = 1,
   lab = L"\Omega_{th}~Data",
)
p = plot!(
   redshift,
   Ωth,
   c = 1,
   ribbon = (Ωth - Ωthl, Ωthu - Ωth),
   lab = "Halo Model",
)
p = plot!(
   twinx(),
   redshift,
   0.2 * Ωth / 1.7755e-8 * 0.049 / Ωb,
   c = 1,
   lw = 2,
   ylab = "Density Weighted Temperature [keV]",
   yaxis = :log,
   legend = false,
   xaxis = false,
   ylims = [1e-9, 2e-6] * 0.2 / 1.7755e-8 * 0.049 / Ωb,
   xlims = [0, 1.5],
   yticks = ([1e-2, 1e-1, 1, 10], ["0.01", "0.1", "1", "10"]),
)
savefig("figure1.pdf")
display(p)

# %% Clean CLASS (the equivalent of the struct_free() in the `main`
# of CLASS. This step is primordial when running in a loop over different
# cosmologies, as you will saturate your memory very fast if you ommit
# it.
cosmo.struct_cleanup()
# If you want to change completely the cosmology, you should also
# clean the arguments, otherwise, if you are simply running on a loop
# of different values for the same parameters, this step is not needed
cosmo.empty()
