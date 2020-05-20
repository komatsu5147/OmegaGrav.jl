function dogravdlnMh(
   pk,
   z::Real,
   Ωm::Real,
   lnMh::Real,
   Ωcb = Ωm;
   virial = false,
   t10MF = false,
)
   ρc = 2.775e11 # in units of h^2 M⊙/Mpc^3
   GN = 3 / 8π / 2998^2 / ρc # = 4.786e-20 Mpc/M⊙ is Newton's constant
   if virial # Determine Δm used by the mass function
      E2 = Ωm * (1 + z)^3 + 1 - Ωm # flat Universe
      Ωmz = Ωm * (1 + z)^3 / E2
      Δvir = 18 * π^2 + 82 * (Ωmz - 1) - 39 * (Ωmz - 1)^2
      Δm = Δvir / Ωmz
   else
      Δm = 200
   end
   spl = dndlnMh(pk, z, Δm, Ωcb, lnMh, t10MF = t10MF)
   # %% Compute ρgrav = (1/2) * \int dlnM dn/dlnM Ag*GM^2/R, in units of h^2 M⊙/Mpc^3
   # c(M) relation is due to Equation (2) and "NFW, F, 0-2" of Table 1 in
   # Duffy et al., MNRAS, 390, L64 (2008)
   if virial # For Δ = Δvir(z)
      RΔh = cbrt(exp(lnMh) * 3 / 4π / (ρc * E2 * Δvir)) # in h^-1 Mpc
      A0, B0, C0 = 7.85, -0.081, -0.71 # Avir, Bvir, Cvir in Table 1 of Duffy et al.
   else # For Δ = Δm = 200
      RΔh = cbrt(exp(lnMh) * 3 / 4π / (ρc * Ωcb * Δm)) / (1 + z) # in h^-1 Mpc
      A0, B0, C0 = 10.14, -0.081, -1.01 # Amean, Bmean, Cmean in Table 1 of Duffy et al.
   end
   c = A0 * (exp(lnMh) / 2e12)^B0 * (1 + z)^C0
   dρdlnMh = 0.5 * spl * GN * exp(2 * lnMh) / RΔh * Ag(c)
   halfW = -dρdlnMh / ρc / Ωcb
   Ωgrav = halfW * Ωcb
end

function dothermdlnMh(
   pk,
   z::Real,
   Ωm::Real,
   h0::Real,
   lnMh::Real,
   Ωcb = Ωm;
   massbias = 1.266,
   t10MF = false,
)
   ρc = 2.775e11 # in units of h^2 M⊙/Mpc^3
   ρceVcm3 = 1.05375e4 # in units of h^2 eV/cm^3
   Δc = 500 # 500 times the critical density of the Universe
   E2 = Ωm * (1 + z)^3 + 1 - Ωm  # flat Universe
   Δm = Δc * E2 / Ωm / (1 + z)^3 # Δm for mass function
   spl = dndlnMh(pk, z, Δm, Ωcb, lnMh, t10MF = t10MF)
   # %% Pressure profile integral
   # Planck's universal pressure profile
   # Reference: Planck Collaboration, A&A, 550, A131 (2013)
   c, γ, α, β, P0 = [1.81, 0.31, 1.33, 4.13, 6.41]
   upp(x) =
      x^2 * P0 * (0.7 / h0)^1.5 * (c * x)^(-γ) * (1 + (c * x)^α)^((γ - β) / α)
   uppint, err = hquadrature(upp, 0, 6)
   # %% Compute ρth = \int dlnM dn/dlnM \int dV Pe, in units of h^2 eV/cm^3
   αp = 0.12
   E2 = Ωm * (1 + z)^3 + 1 - Ωm # for a flat Universe
   RΔh = cbrt(exp(lnMh) * 3 / 4π / (ρc * E2 * Δc)) # in h^-1 Mpc
   Pe =
      1.65 * # eV/cm^3
      (1 / 0.7)^2 * # h^2 cancels h^2 in the critical density in the denominator
      E2^(4 / 3) *
      (exp(lnMh) / massbias / 3e14 / 0.7)^(2 / 3 + αp)
   dρdlnMh = spl * Pe * 4π * RΔh^3 / massbias
   Y = 0.24 # helium abundance
   Ωtherm = (8 - 5Y) / (4 - 2Y) * dρdlnMh * uppint / ρceVcm3
end
