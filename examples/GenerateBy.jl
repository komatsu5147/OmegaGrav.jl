using OmegaGrav
using Dierckx
using Tables, CSV
# %% Compute the matter power spectrum using CLASS
include("compute_pk_class.jl")
# %% Compute Ωtherm at seven redshifts
redshift = [0, 0.3, 0.5, 0.7, 1, 1.3, 1.5]
nred = length(redshift)
b_y = zeros(nred)
for ired = 1:nred
   z = redshift[ired]
   # Define functions to return linear baryon+CDM power spectrum
   # Note: The CLASS code takes wavenumbers in units of 1/Mpc (no h) and
   # return power spectra in units of Mpc^3 (no 1/h^3).
   pkcb_class(kovh) = cosmo.pk_cb_lin(kovh * h0, z) * h0^3
   # Spline interpolate in log(k)
   lnk = log(1e-4):0.05:log(100)
   pkcb = Spline1D(lnk, pkcb_class.(exp.(lnk)))
   # %% Compute the Compton-y-weighted halo bias, b_y
   b_y[ired] = by(x -> pkcb(log(x)), z, Ωm, Ωcb, αp = 0.12)
end

#%% Save table data to by.csv
t = Tables.table([redshift b_y])
header = ["z", "by"]
CSV.write("by.csv", t, header = header)

# %% Clean CLASS (the equivalent of the struct_free() in the `main`
# of CLASS. This step is primordial when running in a loop over different
# cosmologies, as you will saturate your memory very fast if you ommit
# it.
cosmo.struct_cleanup()
# If you want to change completely the cosmology, you should also
# clean the arguments, otherwise, if you are simply running on a loop
# of different values for the same parameters, this step is not needed
cosmo.empty()
