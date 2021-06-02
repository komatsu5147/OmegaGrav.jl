using HaloMF
using HCubature
using Plots
Δm = 200
println("Δm = ", Δm)
A(Δm) = 1 + 0.24log10(Δm) * exp(-(4 / log10(Δm))^4)
a(Δm) = 0.44log10(Δm) - 0.88
C(Δm) = 0.019 + 0.107log10(Δm) + 0.19 * exp(-(4 / log10(Δm))^4)
B, b, c = 0.183, 1.5, 2.4
δc = 1.686
bν(lnν) =
   1 - A(Δm) * exp(0.5 * a(Δm) * lnν) / (exp(0.5 * a(Δm) * lnν) + δc^a(Δm)) +
   B * exp(0.5 * b * lnν) +
   C(Δm) * exp(0.5 * c * lnν)
f0(x) = tinker10MF(x, 0, Δm)
norm1, err = hquadrature(f0, -50, 5)
println("Normalisation(MF) = ", norm1)
f1(x) = tinker10MF(x, 0, Δm) * bν(x)
norm2, err = hquadrature(f1, -50, 5)
println("Normalisation(Bias) = ", norm2)
# %% Plot b(ν)
x = -2:0.05:2
Δm = 200
plot(
   log10.(exp.(x)),
   log10.(bν.(x)),
   xlim = [-0.5, 0.6],
   ylim = [-0.5, 1.2],
   xlabel = "log10(ν)",
   ylabel = "log10(b)",
   ls = :dash,
)
Δm = 1600
plot!(log10.(exp.(x)), log10.(bν.(x)))
