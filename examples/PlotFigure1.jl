using OmegaGrav
using PyCall
using CSV
using Dierckx
using Plots, LaTeXStrings
# %% Call the python wrapper for CLASS, `classy`, via PyCall
classy = pyimport("classy")
# Create an instance of the CLASS wrapper
cosmo = classy.Class()
# Create a dictionary of the cosmological parameters
deg_ncdm = true # 3 massive neutrinos with degenerate mass?
params = Dict(
   "output" => "mPk",
   "P_k_max_h/Mpc" => 100,
   "z_pk" => 5,
   "non linear" => "halofit",
   "A_s" => 2.097e-9,
   "n_s" => 0.9652,
   "k_pivot" => 0.05,
   "h" => 0.6737,
   "omega_b" => 0.02233,
   "omega_cdm" => 0.1198,
   "N_ncdm" => 1,
)
if deg_ncdm
   push!(params, "m_ncdm" => 0.02, "deg_ncdm" => 3, "N_ur" => 0.00641)
   mν = params["m_ncdm"] * params["deg_ncdm"]
else
   push!(params, "m_ncdm" => 0.06, "N_ur" => 2.0328)
   mν = params["m_ncdm"]
end
h0 = params["h"]
Ωcb = (params["omega_b"] + params["omega_cdm"]) / h0^2
Ων = mν / 93.14 / h0^2 #  Ωνh^2 = ∑mν / (93.14 eV)
Ωm = Ωcb + Ων
# Set the parameters and run the CLASS code
cosmo.set(params)
cosmo.compute()

# %% Compute Ωgrav and Ωtherm at redshifts of the data points
d = CSV.read("data/d16_Omega_th_data.csv")
nred = size(d)[1]
redshift = zeros(nred + 1)
Ωghalo = zeros(nred + 1)
Ωgpk = zeros(nred + 1)
Ωth = zeros(nred + 1)
Ωthu = zeros(nred + 1)
Ωthl = zeros(nred + 1)
ΩthB1 = zeros(nred + 1)
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
   lnk = log(1e-4):0.1:log(100)
   pkcb = Spline1D(exp.(lnk), pkcb_class.(exp.(lnk)))
   pknl = Spline1D(exp.(lnk), pknl_class.(exp.(lnk)))
   # %% Compute Ωgrav from non-linear total matter P(k)
   Ωgpk[ired] = ograv_pk(pknl, z, Ωm)
   # %% Compute Ωgrav from Halos, excluding the neutrino contribution
   Ωghalo[ired] = ograv_halo(pkcb, z, Ωm, Ωcb)
   # %% Compute Ωtherm from Halos, excluding the neutrino contribution
   Ωth[ired] = otherm_upp(pkcb, z, Ωm, h0, Ωcb)
   # %% Compute Ωtherm for upper and lower 68% confidence level in B
   Ωthl[ired] = otherm_upp(pkcb, z, Ωm, h0, Ωcb, massbias = 1.315)
   Ωthu[ired] = otherm_upp(pkcb, z, Ωm, h0, Ωcb, massbias = 1.221)
   # %% Compute Ωtherm for no mass bias case for comparison
   ΩthB1[ired] = otherm_upp(pkcb, z, Ωm, h0, Ωcb, massbias = 1)
end

#%% Plot results and save to figure1.pdf
fb = params["omega_b"] / (params["omega_cdm"] + params["omega_b"])
p = scatter(
   (d.z, d.Omega_th),
   yerror = (d.Omega_th .- d.Omega_th_low, d.Omega_th_up .- d.Omega_th),
   ms = 5,
   xlims = [0, 1.5],
   ylims = [1e-9, 1e-6],
   yaxis = :log,
   xlab = "Redshift, z",
   ylab = "Comoving Energy Density Parameters",
   lab = L"\Omega_{th}~Data",
   legend = :topright,
   legendfontsize = 12,
   labelfontsize = 14,
)
p = plot!(
   redshift,
   Ωth,
   ribbon = (Ωth - Ωthl, Ωthu - Ωth),
   lab = L"B=1.266_{-0.045}^{+0.049}",
)
p = plot!(redshift, -Ωgpk, c = :black, lab = L"-\Omega_{grav}^{total}", lw = 2)
p = plot!(redshift, -Ωghalo, c = :green, lab = L"-\Omega_{grav}^{halo}", lw = 2)
p = plot!(
   redshift,
   -fb * Ωghalo,
   c = :green,
   lab = L"-f_b\Omega_{grav}^{halo}",
   ls = :dashdot,
   lw = 2,
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
