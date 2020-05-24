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
   # Spline interpolate in log(k)
   lnk = log(1e-4):0.05:log(100)
   pkcb = Spline1D(lnk, pkcb_class.(exp.(lnk)))
   # %% Compute Ωgrav from Halos, excluding the neutrino contribution
   Ωghalo[ired] = ograv_halo(x -> pkcb(log(x)), z, Ωm, Ωcb)
   # %% Compute Ωtherm from Halos, excluding the neutrino contribution
   Ωth[ired] = otherm_upp(x -> pkcb(log(x)), z, Ωm, h0, Ωcb)
   # %% Compute Ωtherm for upper and lower 68% confidence level in B
   Ωthl[ired] = otherm_upp(x -> pkcb(log(x)), z, Ωm, h0, Ωcb, massbias = 1.315)
   Ωthu[ired] = otherm_upp(x -> pkcb(log(x)), z, Ωm, h0, Ωcb, massbias = 1.221)
   # %% Compute Ωtherm for no mass bias case (B=1) for comparison
   ΩthB1[ired] = otherm_upp(x -> pkcb(log(x)), z, Ωm, h0, Ωcb, massbias = 1)
   # %% Compute Ωtherm for Komatsu-Seljak profile
   # ΩthB1[ired] = otherm_ks(x -> pkcb(log(x)), z, Ωm, 0.157, Ωcb)
end

#%% Plot results and save to figure3.pdf
fb = params["omega_b"] / (params["omega_cdm"] + params["omega_b"])
ii = findall(x -> x < 1, d.z)
p = scatter(
   d.z[ii],
   -d.Omega_th[ii] ./ Ωghalo[ii.+1] / fb,
   yerror = (
      -(d.Omega_th .- d.Omega_th_low) ./ Ωghalo[2:end] / fb,
      -(d.Omega_th_up .- d.Omega_th) ./ Ωghalo[2:end] / fb,
   ),
   ms = 5,
   c = 1,
   #xticks = (1:4, string.(1:4)),
   ylims = [0.1, 1.6],
   xlab = "Redshift, z",
   ylab = L"\Omega_{th}/(-f_b\Omega_{grav}^{halo})",
   lab = "Data",
   legend = :topright,
   legendfontsize = 12,
   labelfontsize = 15,
)
ii = findall(x -> x > 1, d.z)
u = zeros(length(ii))
v = -0.1 * ones(length(ii))
p = quiver!(
   d.z[ii],
   -d.Omega_th_up[ii] ./ Ωghalo[ii.+1] / fb,
   quiver = (u, v),
   c = 1,
)
p = plot!(
   redshift,
   -Ωth ./ Ωghalo / fb,
   ribbon = (-(Ωth - Ωthl) ./ Ωghalo / fb, -(Ωthu - Ωth) ./ Ωghalo / fb),
   lab = L"B=1.27_{-0.04}^{+0.05}",
   lw = 5,
   c = 1,
)
p = plot!(redshift, -ΩthB1 ./ Ωghalo / fb, ls = :dash, lab = L"B=1", lw = 5)
savefig("figure3.pdf")
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
