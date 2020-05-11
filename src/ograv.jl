"""
   ograv_pk(pk, z, Ωm; kmin = 5e-4, kmax = 1e2)

Comoving density parameter of gravitational binding energy computed from a matter power spectrum.

*Reference*: Equation (60) of Fukugita & Peebles, ApJ, 616, 643 (2004)
- With the redshift evolution factor of Chiang, Makiya, Komatsu & Ménard (in prep)

# Arguments
- `pk`(k): a function which returns a matter power spectrum with the argument k being the comoving wavenumber.
    - **pk times k^3 must be dimensionless**. For example, if k is in units of h/Mpc, `pk` must be in units of Mpc^3/h^3.
- `z::Real`: redshift.
- `Ωm::Real`: present-day matter density parameter.

# Optional arguments
- `kmin::Real=5e-4`: minimum wavenumber for integration, ``∫_{kmin}^{kmax} dk P(k)``.
- `kmax::Real=1e2`: maximum wavenumber for integration, ``∫_{kmin}^{kmax} dk P(k)``.
"""
function ograv_pk(pk, z::Real, Ωm::Real; kmin = 5e-4, kmax = 1e2)
   res, err = hquadrature(pk, kmin, kmax)
   halfW = -3 * Ωm * (1 + z) / 2998^2 / 16 / π^2 * res
   Ωgrav = halfW * Ωm
end

"""
   ograv_halo(pk, z, Ωm; Mmin = 5e8, Mmax = 5e15, virial = false, t10MF = false)

Comoving density parameter of gravitational binding energy computed from a matter power spectrum.

*Reference*: Equation (TBD) of Chiang, Makiya, Komatsu & Ménard (in prep)

# Arguments
- `pk`(k): a function which returns a matter power spectrum with the argument k being the comoving wavenumber.
    - **pk times k^3 must be dimensionless**. For example, if k is in units of h/Mpc, `pk` must be in units of Mpc^3/h^3.
- `z::Real`: redshift.
- `Ωm::Real`: present-day matter density parameter.

# Optional arguments
- `Mmin::Real=5e8`: minimum mass for integration, ``∫_{Mmin}^{Mmax} dM dn/dM Ag GM^2/R``.
- `Mmax::Real=5e15`: maximum mass for integration, ``∫_{Mmin}^{Mmax} dM dn/dM Ag GM^2/R``.
- `virial::Bool=false`: if `true`, use the virial overdensity `Δvir`. If `false` (the default), use `Δm=200`.
- `t10MF::Bool=false`: if `true`, use `tinker10MF` for the halo multiplicity function. If `false` (the default), use `tinker08MF`.
"""
function ograv_halo(
   pk,
   z::Real,
   Ωm::Real;
   Mmin = 5e8,
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
   lnMh = range(log(Mmin), length = nmass, log(Mmax))
   dndlnMh = zeros(nmass)
   for i = 1:nmass
      Rh = cbrt(exp(lnMh[i]) * 3 / 4π / ρc / Ωm) # in units of Mpc/h
      σ2 = sigma2(pk, Rh)
      dlnσ2dlnRh = Rh * dsigma2dR(pk, Rh) / σ2
      lnν = 2 * log(1.6865) - log(σ2)
      if t10MF
         MF = tinker10MF(lnν, z, Δm)
      else
         MF = tinker08MF(lnν, z, Δm)
      end
      dndlnMh[i] = -dlnσ2dlnRh * MF / 4π / Rh^3 # in units of h^3 Mpc^-3
   end
   spl = Spline1D(lnMh, dndlnMh)
   # %% Compute ρgrav = (1/2) * \int dlnM dn/dlnM Ag*GM^2/R, in units of h^2 M⊙/Mpc^3
   # c(M) relation is due to Equation (2) and "NFW, F, 0-2" of Table 1 in
   # Duffy et al., MNRAS, 390, L64 (2008)
   function dρdlnMh(lnMh)
      if virial # For Δ = Δvir(z)
         RΔh = cbrt(exp(lnMh) * 3 / 4π / (ρc * E2 * Δvir)) # in h^-1 Mpc
         A, B, C = 7.85, -0.081, -0.71 # Avir, Bvir, Cvir in Table 1 of Duffy et al.
      else # For Δ = Δm = 200
         RΔh = cbrt(exp(lnMh) * 3 / 4π / (ρc * Ωm * Δm)) / (1 + z) # in h^-1 Mpc
         A, B, C = 10.14, -0.081, -1.01 # Amean, Bmean, Cmean in Table 1 of Duffy et al.
      end
      c = A * (exp(lnMh) / 2e12)^B / (1 + z)^1.01
      dρdlnMh = 0.5 * spl(lnMh) * GN * exp(2 * lnMh) / RΔh * Ag(c)
   end
   ρgrav, err = hquadrature(dρdlnMh, log(Mmin), log(Mmax))
   halfW = -ρgrav / ρc / Ωm
   Ωgrav = halfW * Ωm
end

"""
   Ag(c:Real; method = "poly4")

Integral of the squared Fourier transform of an NFW density profile, ``(1/π)∫_0^∞ dx |u(x,c)|^2``.

*Reference*: Equation (TBD) of Chiang, Makiya, Komatsu & Ménard (in prep)

# Arguments
- `c::Real`: concentration parameter.

# Optional arguments
- `method::String="poly4"`: approximation method.
   - `poly4` (the default) or other values for 4th order polynomial fit.
   - `poly3` for 3rd order polynomial fit.
   - `exact` for numerical integration.
"""
function Ag(c::Real; method = "poly4")
   if method == "poly3" # 3rd order polynomial fit
      y = c - 5
      Ag = 1.0184 + 5.389e-2 * y - 1.514e-3 * y^2 + 3.845e-5 * y^3
   elseif method == "exact" # exact integration
      f = log(1 + c) - c / (1 + c)
      u2(x) =
         (
            sin(x / c) * (sinint((1 / c + 1) * x) - sinint(x / c)) +
            cos(x / c) * (cosint((1 / c + 1) * x) - cosint(x / c)) -
            sin(x) / (1 / c + 1) / x
         )^2
      res, err = hquadrature(u2, 0, 100)
      Ag = res / π / f^2
   else # 4th order polynomial fit
      y = c - 5
      Ag =
         1.0202 + 5.376e-2 * y - 1.842e-3 * y^2 + 9.489e-5 * y^3 -
         2.352e-6 * y^4
   end
end