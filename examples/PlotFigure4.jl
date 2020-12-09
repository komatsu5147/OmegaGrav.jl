using OmegaGrav
using MatterPower
using Dierckx
using Plots, LaTeXStrings
using OrdinaryDiffEq
# %% Compute the matter power spectrum using CLASS
include("compute_pk_class.jl")
# %% Compute the analytical linear theory result for K/W
# Reference: Equation (12) of Davis, Miller, White, ApJ, 490, 63 (1997).
zini = 9
a = 1/(1+zini):0.01:1.0 # scale factor
nred = length(a)
sol_δθ = setup_growth(Ωm, 1 - Ωm)
KovW = zeros(nred)
for i = 1:nred
    KovW[i] = -2 / 3 * (sol_δθ(a[i])[2] / sol_δθ(a[i])[1])^2 * a[i] / Ωm
end

# %% Compute `Ωgrav = Ωm * W / 2` at various redshifts
Ωgpk = zeros(nred)
Ωglin = zeros(nred)
for ired = 1:nred
    z = 1 / a[ired] - 1
    # Define functions to return linear and non-linear power spectra
    # Note: The CLASS code takes wavenumbers in units of 1/Mpc (no h) and
    # return power spectra in units of Mpc^3 (no 1/h^3).
    pknl_class(kovh) = cosmo.pk(kovh * h0, z) * h0^3
    pklin_class(kovh) = cosmo.pk_lin(kovh * h0, z) * h0^3
    # Spline interpolate in log(k)
    lnk = log(1e-4):0.05:log(100)
    pknl = Spline1D(lnk, pknl_class.(exp.(lnk)))
    pklin = Spline1D(lnk, pklin_class.(exp.(lnk)))
    # %% Compute Ωgrav from non-linear total matter P(k)
    Ωgpk[ired] = ograv_pk(x -> pknl(log(x)), z, Ωm)
    # %% Compute Ωgrav from linear total matter P(k)
    Ωglin[ired] = ograv_pk(x -> pklin(log(x)), z, Ωm)
end
Ωgrav_pklin = Spline1D(a, Ωglin)
Ωgrav_pknl = Spline1D(a, Ωgpk)

# %% Solve the Layzer-Irvine equation
tspan = (1 / (1 + zini), 1.0)
t0 = tspan[1]
## Non-linear P(k)
W(x) = 2 * Ωgrav_pknl(x)
dWda(x) = 2 * derivative(Ωgrav_pknl, x)
f(u, p, t) = -(2 * u + W(t)) / t - dWda(t)
u0 = KovW[1] * W(t0)
prob = ODEProblem(f, u0, tspan)
sol_pknl = solve(prob, Tsit5())
## Linear P(k): This result must agree with the analytical linear theory result
W(x) = 2 * Ωgrav_pklin(x)
dWda(x) = 2 * derivative(Ωgrav_pklin, x)
f(u, p, t) = -(2 * u + W(t)) / t - dWda(t)
u0 = KovW[1] * W(t0)
prob = ODEProblem(f, u0, tspan)
sol_pklin = solve(prob, Tsit5())

#%% Plot results and save to figure4.pdf
## Plot against the scale factor
#p = plot(
#    a,
#    sol_pknl.(a)[:, 1],
#    lw = 2,
#    xlab = "Scale factor, a",
#    ylab = "Comiving Density Parameters",
#    lab = L"\Omega_K: Nonlinear",
#    legend = :topleft,
#    legendfontsize = 12,
#    labelfontsize = 15,
#    ylims = [1e-7, 8e-7],
#    xlims = [0.2, 1.0],
#)
#p = plot!(a, -Ωgrav_pknl.(a), lw = 2, lab = L"-\Omega_W/2: Nonlinear")
#p = plot!(a, sol_pklin.(a)[:, 1], c = 1, lw = 2, ls = :dot, lab = L"\Omega_K: Linear")
#p = plot!(a, -Ωgrav_pklin.(a), c = 2, lw = 2, ls = :dot, lab = L"-\Omega_W/2: Linear")
#p = plot!(a, KovW .* W.(a), ls = :dot, lab = L"K: Linear,~Analytical")
## Plot against the redshift
redshift = 1 ./ a .- 1
p = plot(
    redshift,
    sol_pknl.(a)[:, 1],
    lw = 2,
    xlab = "Redshift, z",
    ylab = "Comiving Density Parameters",
    lab = L"\Omega_K: Nonlinear",
    legend = :topright,
    legendfontsize = 12,
    labelfontsize = 15,
    ylims = [1e-7, 8e-7],
    xlims = [0, 3.0],
)
p = plot!(redshift, -Ωgrav_pknl.(a), lw = 2, lab = L"-\Omega_W/2: Nonlinear")
p = plot!(
    redshift,
    sol_pklin.(a)[:, 1],
    c = 1,
    lw = 2,
    ls = :dot,
    lab = L"\Omega_K: Linear",
)
p = plot!(
    redshift,
    -Ωgrav_pklin.(a),
    c = 2,
    lw = 2,
    ls = :dot,
    lab = L"-\Omega_W/2: Linear",
)
savefig("figure4.pdf")
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
