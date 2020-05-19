"""
   otherm_upp(pk, z, Ωm, h0[, Ωcb=Ωm]; massbias = 1.266, Mmin = 1e11, Mmax = 5e15, t10MF = false)

Comoving density parameter of thermal energy of ionized baryonsgas in collapsed structures using Planck's universal pressure profile.

*Reference*: Equation (TBD) of Chiang, Makiya, Komatsu & Ménard (in prep)

# Arguments
- `pk`(k): a function which returns a matter power spectrum with the argument k being the comoving wavenumber.
   - `pk` is in units of (Mpc/h)^3 and the argument k is in units of h/Mpc.
- `z::Real`: redshift.
- `Ωm::Real`: present-day total matter density parameter.
- `h0::Real`: present-day Hubble constant in units of 100 km/s/Mpc.

# Optional Arguments
- `Ωcb::Real=Ωm`: present-day baryon + cold dark matter density parameter.

# Optional keyword arguments
- `massbias::Real=1.266`: mass bias parameter for Planck's universal pressure profile, `B=Mtrue/Mplanck`.
- `Mmin::Real=1e11`: minimum mass for integration, ``∫_{Mmin}^{Mmax} dn/dM ∫dV Pe``.
- `Mmax::Real=5e15`: maximum mass for integration, ``∫_{Mmin}^{Mmax} dn/dM ∫dV Pe``.
- `t10MF::Bool=false`: if `true`, use `tinker10MF` for the halo multiplicity function. If `false` (the default), use `tinker08MF`.
"""
function otherm_upp(
   pk,
   z::Real,
   Ωm::Real,
   h0::Real,
   Ωcb = Ωm;
   massbias = 1.266,
   Mmin = 1e11,
   Mmax = 5e15,
   t10MF = false,
)
   ρc = 2.775e11 # in units of h^2 M⊙/Mpc^3
   ρceVcm3 = 1.05375e4 # in units of h^2 eV/cm^3
   Δc = 500 # 500 times the critical density of the Universe
   E2 = Ωm * (1 + z)^3 + 1 - Ωm  # flat Universe
   Δm = Δc * E2 / Ωm / (1 + z)^3 # Δm for mass function
   nmass = 100
   func = zeros(nmass)
   lnMh = range(log(Mmin), length = nmass, log(Mmax))
   for i = 1:nmass
      func[i] = dndlnMh(pk, z, Δm, Ωcb, lnMh[i], t10MF = t10MF)
   end
   spl = Spline1D(lnMh, func)
   # %% Pressure profile integral
   # Planck's universal pressure profile
   # Reference: Planck Collaboration, A&A, 550, A131 (2013)
   c, γ, α, β, P0 = [1.81, 0.31, 1.33, 4.13, 6.41]
   upp(x) =
      x^2 * P0 * (0.7 / h0)^1.5 * (c * x)^(-γ) * (1 + (c * x)^α)^((γ - β) / α)
   uppint, err = hquadrature(upp, 0, 6)
   # %% Compute ρth = \int dlnM dn/dlnM \int dV Pe, in units of h^2 eV/cm^3
   function dρdlnMh(lnMh)
      αp = 0.12
      E2 = Ωm * (1 + z)^3 + 1 - Ωm # for a flat Universe
      RΔh = cbrt(exp(lnMh) * 3 / 4π / (ρc * E2 * Δc)) # in h^-1 Mpc
      Pe =
         1.65 * # eV/cm^3
         (1 / 0.7)^2 * # h^2 cancels h^2 in the critical density in the denominator
         E2^(4 / 3) *
         (exp(lnMh) / massbias / 3e14 / 0.7)^(2 / 3 + αp)
      dρdlnMh = spl(lnMh) * Pe * 4π * RΔh^3 / massbias
   end
   ρe, err = hquadrature(dρdlnMh, log(Mmin), log(Mmax))
   Y = 0.24 # helium abundance
   Ωtherm = (8 - 5Y) / (4 - 2Y) * ρe * uppint / ρceVcm3
end

"""
   otherm_ks(pk, z, Ωm, fb[, Ωcb=Ωm]; Mmin = 1e11, Mmax = 5e15, virial = false, t10MF = false)

Comoving density parameter of thermal energy of ionized baryonsgas in collapsed structures using Komatsu-Seljak pressure profile.

*Reference*: Equation (TBD) of Chiang, Makiya, Komatsu & Ménard (in prep)

# Arguments
- `pk`(k): a function which returns a matter power spectrum with the argument k being the comoving wavenumber.
    - **pk times k^3 must be dimensionless**. For example, if k is in units of h/Mpc, `pk` must be in units of Mpc^3/h^3.
- `z::Real`: redshift.
- `Ωm::Real`: present-day total matter density parameter.
- `fb::Real`: mean baryon fraction.

# Optional Arguments
- `Ωcb::Real=Ωm`: present-day baryon + cold dark matter density parameter.

# Optional keyword arguments
- `Mmin::Real=1e11`: minimum mass for integration, ``∫_{Mmin}^{Mmax} dn/dM ∫dV Pe``.
- `Mmax::Real=5e15`: maximum mass for integration, ``∫_{Mmin}^{Mmax} dn/dM ∫dV Pe``.
- `virial::Bool=false`: if `true`, use the virial overdensity `Δvir`. If `false` (the default), use `Δm=200`.
- `t10MF::Bool=false`: if `true`, use `tinker10MF` for the halo multiplicity function. If `false` (the default), use `tinker08MF`.
"""
function otherm_ks(
   pk,
   z::Real,
   Ωm::Real,
   fb::Real,
   Ωcb = Ωm;
   Mmin = 1e11,
   Mmax = 5e15,
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
   nmass = 100
   func = zeros(nmass)
   lnMh = range(log(Mmin), length = nmass, log(Mmax))
   for i = 1:nmass
      func[i] = dndlnMh(pk, z, Δm, Ωcb, lnMh[i], t10MF = t10MF)
   end
   spl = Spline1D(lnMh, func)
   ## Pressure profile integral
   # Komatsu-Seljak Pressure Profile
   # Reference: Komatsu & Seljak, MNRAS, 327, 1353 (2001)
   # Fitting formulae for γ(c) and η(c) are valid for 1<c<25.
   # Reference: Equations (17,18) of Komatsu & Seljak, MNRAS, 336, 1256 (2002)
   γ(c) = 1.137 + 8.94e-2 * log(c / 5) - 3.68e-3 * (c - 5)
   η0(c) = 2.235 + 0.202 * (c - 5) - 1.16e-3 * (c - 5)^2
   B(c) = 3 / η0(c) * (γ(c) - 1) / γ(c) / (log(1 + c) / c - 1 / (1 + c))
   ygas(x, c) = (1 - B(c) * (1 - log(1 + x) / x))^(1 / (γ(c) - 1))
   ksint = zeros(25)
   for c = 1:25
      ks(x) = x^2 * ygas(x, c)^γ(c)
      ksint[c], err = hquadrature(ks, 0, 3c) # integrated out to r = 3*r200m
   end
   c = 1:25
   ksints = Spline1D(c, ksint)
   # %% Compute ρth = \int dlnM dn/dlnM \int dV Pe, in units of h^2 eV/cm^3
   ## Gas thermal energy density: Komatsu-Seljak profile
   function dρdlnMh(lnMh)
      if virial # For Δ = Δvir(z)
         RΔh = cbrt(exp(lnMh) * 3 / 4π / (ρc * E2 * Δvir)) # in h^-1 Mpc
         A0, B0, C0 = 7.85, -0.081, -0.71 # Avir, Bvir, Cvir in Table 1 of Duffy et al.
      else # For Δ = Δm = 200
         RΔh = cbrt(exp(lnMh) * 3 / 4π / (ρc * Ωcb * Δm)) / (1 + z) # in h^-1 Mpc
         A0, B0, C0 = 10.14, -0.081, -1.01 # Amean, Bmean, Cmean in Table 1 of Duffy et al.
      end
      c = A0 * (exp(lnMh) / 2e12)^B0 * (1 + z)^C0
      mc = log(1 + c) - c / (1 + c)
      ρgnorm = fb * (c * (1 + c)^2 * mc * ygas(c, c))^-1 # = 4πρgas(0)rs^3/Mgas
      dρdlnMh =
         spl(lnMh) * GN * exp(2 * lnMh) / RΔh * η0(c) / 3 * ρgnorm * ksints(c) # in units of h^2 M⊙/Mpc^3
   end
   res, err = hquadrature(dρdlnMh, log(Mmin), log(Mmax))
   Ωtherm = res / ρc
end
