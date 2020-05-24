using OmegaGrav
using PyCall
using Dierckx
using Plots, LaTeXStrings
include("dodlnMh.jl")
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

# %% Compute dΩgrav/dlnMh and dΩtherm/dlnMh at z = 0
z = 0
# Define functions to return linear and non-linear power spectra
# Note: The CLASS code takes wavenumbers in units of 1/Mpc (no h) and
# return power spectra in units of Mpc^3 (no 1/h^3).
pkcb_class(kovh) = cosmo.pk_cb_lin(kovh * h0, z) * h0^3
# Spline interpolate in log(k)
lnk = log(1e-4):0.1:log(100)
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
   -fb * dΩgdlnMh.(lnMh),
   xaxis = :log,
   yaxis = :log,
   xlab = "Halo Mass [M⊙/h]",
   ylab = L"d\Omega/d\ln M",
   lab = L"-f_b\Omega_{grav}^{halo}",
   legend = :topleft,
   legendfontsize = 12,
   labelfontsize = 15,
   lw = 5,
   c = :green,
   ls = :dashdot,
   ylims = [1e-10, 1e-8],
   xlims = [1e11, 5e15],
)
p = plot!(exp.(lnMh), dΩthdlnMh.(lnMh), lw = 5, c = 1, lab = L"\Omega_{th}")

# %% Compute dΩgrav/dlnMh and dΩtherm/dlnMh at z = 1
z = 1
pkcb_class(kovh) = cosmo.pk_cb_lin(kovh * h0, z) * h0^3
pkcb = Spline1D(lnk, pkcb_class.(exp.(lnk)))
dΩgdlnMh(lnMh) = -fb * dogravdlnMh(x -> pkcb(log(x)), z, Ωm, lnMh, Ωcb)
dΩthdlnMh(lnMh) = dothermdlnMh(x -> pkcb(log(x)), z, Ωm, h0, lnMh, Ωcb)
p = plot!(
   exp.(lnMh),
   dΩgdlnMh.(lnMh),
   lw = 2,
   c = :green,
   ls = :dashdot,
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
