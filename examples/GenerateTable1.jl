using OmegaGrav
using Dierckx
using Tables, CSV
# %% Compute the matter power spectrum using CLASS
include("compute_pk_class.jl")
# %% Compute Ωgrav at seven redshifts
redshift = [0, 0.3, 0.5, 0.7, 1, 1.3, 1.5]
nred = length(redshift)
Ωgpk = zeros(nred)
Ωglin = zeros(nred)
Ωghalo = zeros(nred)
Ωgh12 = zeros(nred)
Ωgh13 = zeros(nred)
Ωgh14 = zeros(nred)
# ΩghaloT10 = zeros(nred)
# Ωgh12T10 = zeros(nred)
# Ωgh13T10 = zeros(nred)
# Ωgh14T10 = zeros(nred)
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
   # %% Compute Ωgrav from Halos, excluding the neutrino contribution, with Tinker10MF
   # ΩghaloT10[ired] = ograv_halo(x -> pkcb(log(x)), z, Ωm, Ωcb, t10MF = true)
   # Ωgh12T10[ired] = ograv_halo(x -> pkcb(log(x)), z, Ωm, Ωcb, Mmin = 1e12, t10MF = true)
   # Ωgh13T10[ired] = ograv_halo(x -> pkcb(log(x)), z, Ωm, Ωcb, Mmin = 1e13, t10MF = true)
   # Ωgh14T10[ired] = ograv_halo(x -> pkcb(log(x)), z, Ωm, Ωcb, Mmin = 1e14, t10MF = true)
end

#%% Save table data to table1.csv
fb = params["omega_b"] / (params["omega_cdm"] + params["omega_b"])
t = Tables.table([redshift Ωglin Ωgpk Ωghalo Ωgh12 Ωgh13 Ωgh14 2 * fb * Ωghalo / 3 2 * fb * (Ωgpk - Ωglin) / 3])
header = [
   "z",
   "Omega_W/2_pklin",
   "Omega_W/2_pknl",
   "Omega_W/2_halo",
   "Omega_W/2_halo_Mgt1e12",
   "Omega_W/2_halo_Mgt1e13",
   "Omega_W/2_halo_Mgt1e14",
   "Omega_W/2_halo_times_twothirdfb",
   "Omega_W/2_pkdiff_times_twothirdfb",
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
