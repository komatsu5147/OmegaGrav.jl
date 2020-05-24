using OmegaGrav
using PyCall
using Dierckx
using Tables, CSV
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

# %% Compute Ωgrav and Ωtherm at seven redshifts
redshift = [0, 0.3, 0.5, 0.7, 1, 1.3, 1.5]
nred = length(redshift)
Ωgpk = zeros(nred)
Ωglin = zeros(nred)
Ωghalo = zeros(nred)
Ωgh12 = zeros(nred)
Ωgh13 = zeros(nred)
Ωgh14 = zeros(nred)
for ired = 1:nred
   z = redshift[ired]
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
   Ωgh12[ired] = ograv_halo(x -> pkcb(log(x)), z, Ωm, Ωcb, Mmin = 1e12)
   Ωgh13[ired] = ograv_halo(x -> pkcb(log(x)), z, Ωm, Ωcb, Mmin = 1e13)
   Ωgh14[ired] = ograv_halo(x -> pkcb(log(x)), z, Ωm, Ωcb, Mmin = 1e14)
end

#%% Save table data to table1.csv
fb = params["omega_b"] / (params["omega_cdm"] + params["omega_b"])
t = Tables.table([redshift Ωglin Ωgpk Ωghalo Ωgh12 Ωgh13 Ωgh14 fb * Ωghalo])
header = [
   "z",
   "Omega_grav_pklin",
   "Omega_grav_pknl",
   "Omega_grav_halo",
   "Omega_grav_halo_Mgt1e12",
   "Omega_grav_halo_Mgt1e13",
   "Omega_grav_halo_Mgt1e14",
   "Omega_grav_halo_times_fb",
]
CSV.write("table1.csv", t, header = header)

# %% Clean CLASS (the equivalent of the struct_free() in the `main`
# of CLASS. This step is primordial when running in a loop over different
# cosmologies, as you will saturate your memory very fast if you ommit
# it.
cosmo.struct_cleanup()
# If you want to change completely the cosmology, you should also
# clean the arguments, otherwise, if you are simply running on a loop
# of different values for the same parameters, this step is not needed
cosmo.empty()
