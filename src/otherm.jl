"""
   otherm_upp(pk, z, Ωm, h0[, Ωcb=Ωm]; massbias = 1.266, Mmin = 1e11, Mmax = 5e15, t10MF = false)

Comoving density parameter of thermal energy of ionized baryonsgas in collapsed structures using Planck's universal pressure profile

*Reference*: Equation (TBD) of Chiang, Makiya, Komatsu & Ménard (in prep)

# Arguments
- `pk`(k): a function which returns a matter power spectrum with the argument k being the comoving wavenumber.
    - **pk times k^3 must be dimensionless**. For example, if k is in units of h/Mpc, `pk` must be in units of Mpc^3/h^3.
- `z::Real`: redshift.
- `Ωm::Real`: present-day total matter density parameter.
- `h0::Real`: present-day Hubble constant in units of 100 km/s/Mpc.

# Optional Arguments
- `Ωcb::Real=Ωm`: present-day baryon + cold dark matter density parameter.

# Optional keyword arguments
- `massbias::Real=1.266`: mass bias parameter for Planck's universal pressure profile, `B=Mtrue/Mplanck`.
- `Mmin::Real=1e11`: minimum mass for integration, ``∫_{Mmin}^{Mmax} dM dn/dM Ag GM^2/R``.
- `Mmax::Real=5e15`: maximum mass for integration, ``∫_{Mmin}^{Mmax} dM dn/dM Ag GM^2/R``.
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
   nmass = 100
   lnMh = range(log(Mmin), length = nmass, log(Mmax))
   dndlnMh = zeros(nmass)
   for i = 1:nmass
      Rh = cbrt(exp(lnMh[i]) * 3 / 4π / ρc / Ωcb) # in units of Mpc/h
      σ2 = sigma2(pk, Rh)
      dlnσ2dlnRh = Rh * dsigma2dR(pk, Rh) / σ2
      lnν = 2 * log(1.6865) - log(σ2)
      E2 = Ωm * (1 + z)^3 + 1 - Ωm  # flat Universe
      Δm = Δc * E2 / Ωm / (1 + z)^3 # Δm for mass function
      if t10MF
         MF = tinker10MF(lnν, z, Δm)
      else
         MF = tinker08MF(lnν, z, Δm)
      end
      dndlnMh[i] = -dlnσ2dlnRh * MF / 4π / Rh^3 # in units of h^3 Mpc^-3
   end
   spl = Spline1D(lnMh, dndlnMh)
   # %% Pressure profile integral
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
      spl(lnMh) * Pe * 4π * RΔh^3 / massbias
   end
    ρe, err = hquadrature(dρdlnMh, log(Mmin), log(Mmax))
    Y = 0.24 # helium abundance
   Ωtherm = (8-5Y)/(4-2Y) * ρe * uppint / ρceVcm3
end
