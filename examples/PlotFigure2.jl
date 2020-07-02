using OmegaGrav
using Dierckx
using Plots, LaTeXStrings
include("dodlnMh.jl")
# %% Compute the matter power spectrum using CLASS
include("compute_pk_class.jl")
# %% Compute dΩgrav/dlnMh and dΩtherm/dlnMh at z = 0
z = 0
# Define functions to return linear and non-linear power spectra
# Note: The CLASS code takes wavenumbers in units of 1/Mpc (no h) and
# return power spectra in units of Mpc^3 (no 1/h^3).
pkcb_class(kovh) = cosmo.pk_cb_lin(kovh * h0, z) * h0^3
# Spline interpolate in log(k)
lnk = log(1e-4):0.05:log(100)
pkcb = Spline1D(lnk, pkcb_class.(exp.(lnk)))
# %% Define a function to return dΩgrav/dlnMh from Halos, excluding the neutrino contribution
dΩgdlnMh(lnMh) = dogravdlnMh(x -> pkcb(log(x)), z, Ωm, lnMh, Ωcb)
# %% Define a function to return dΩtherm/dlnMh from Halos, excluding the neutrino contribution
dΩthdlnMh(lnMh) = dothermdlnMh(x -> pkcb(log(x)), z, Ωm, h0, lnMh, Ωcb)

#%% Plot results and save to figure2.pdf
fb = params["omega_b"] / (params["omega_cdm"] + params["omega_b"])
lnMh = log(1e11):0.1:log(5e15)
p = plot(
   exp.(lnMh),
   -2 * fb * dΩgdlnMh.(lnMh) / 3,
   xaxis = :log,
   yaxis = :log,
   xlab = "Halo Mass [M⊙/h]",
   ylab = L"d\Omega/d\ln M",
   lab = L"\Omega_W^{ref,halo}",
   legend = :topleft,
   legendfontsize = 12,
   labelfontsize = 15,
   lw = 5,
   c = :green,
   ls = :dot,
   ylims = [1e-10, 1e-8],
   xlims = [1e11, 5e15],
)
p = plot!(exp.(lnMh), dΩthdlnMh.(lnMh), lw = 5, c = 1, lab = L"\Omega_{th}")

# %% Compute dΩgrav/dlnMh and dΩtherm/dlnMh at z = 1
z = 1
pkcb_class(kovh) = cosmo.pk_cb_lin(kovh * h0, z) * h0^3
pkcb = Spline1D(lnk, pkcb_class.(exp.(lnk)))
dΩgdlnMh(lnMh) = dogravdlnMh(x -> pkcb(log(x)), z, Ωm, lnMh, Ωcb)
dΩthdlnMh(lnMh) = dothermdlnMh(x -> pkcb(log(x)), z, Ωm, h0, lnMh, Ωcb)
p = plot!(
   exp.(lnMh),
   -2 * fb * dΩgdlnMh.(lnMh) / 3,
   lw = 2,
   c = :green,
   ls = :dot,
   lab = "",
)
p = plot!(exp.(lnMh), dΩthdlnMh.(lnMh), lw = 2, c = 1, lab = "")
savefig("figure2.pdf")
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
