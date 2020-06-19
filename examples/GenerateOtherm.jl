using OmegaGrav
using Dierckx
using Tables, CSV
# %% Compute the matter power spectrum using CLASS
include("compute_pk_class.jl")
# %% Compute Ωtherm at seven redshifts
redshift = [0, 0.3, 0.5, 0.7, 1, 1.3, 1.5]
nred = length(redshift)
Ωth = zeros(nred)
Ωthu = zeros(nred)
Ωthl = zeros(nred)
ΩthB1 = zeros(nred)
for ired = 1:nred
   z = redshift[ired]
   # Define functions to return linear baryon+CDM power spectrum
   # Note: The CLASS code takes wavenumbers in units of 1/Mpc (no h) and
   # return power spectra in units of Mpc^3 (no 1/h^3).
   pkcb_class(kovh) = cosmo.pk_cb_lin(kovh * h0, z) * h0^3
   # Spline interpolate in log(k)
   lnk = log(1e-4):0.05:log(100)
   pkcb = Spline1D(lnk, pkcb_class.(exp.(lnk)))
   # %% Compute Ωtherm from Halos, excluding the neutrino contribution
   Ωth[ired] = otherm_upp(x -> pkcb(log(x)), z, Ωm, h0, Ωcb)
   # %% Compute Ωtherm for upper and lower 68% confidence level in B
   Ωthl[ired] = otherm_upp(x -> pkcb(log(x)), z, Ωm, h0, Ωcb, massbias = 1.315)
   Ωthu[ired] = otherm_upp(x -> pkcb(log(x)), z, Ωm, h0, Ωcb, massbias = 1.221)
   # %% Compute Ωtherm for no mass bias case for comparison
   ΩthB1[ired] = otherm_upp(x -> pkcb(log(x)), z, Ωm, h0, Ωcb, massbias = 1)
end

#%% Save table data to table1.csv
fb = params["omega_b"] / (params["omega_cdm"] + params["omega_b"])
t = Tables.table([redshift Ωth Ωth - Ωthl Ωthu - Ωth ΩthB1])
header = [
   "z",
   "Omega_th_B1p266",
   "Omega_th_err_low",
   "Omega_th_err_high",
   "Omega_th_B1",
]
CSV.write("otherm_upp.csv", t, header = header)

# %% Clean CLASS (the equivalent of the struct_free() in the `main`
# of CLASS. This step is primordial when running in a loop over different
# cosmologies, as you will saturate your memory very fast if you ommit
# it.
cosmo.struct_cleanup()
# If you want to change completely the cosmology, you should also
# clean the arguments, otherwise, if you are simply running on a loop
# of different values for the same parameters, this step is not needed
cosmo.empty()
