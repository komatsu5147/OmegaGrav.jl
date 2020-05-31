"""
    onehalo(pk, z, Ωm, k[, Ωcb]; Mmin = 1e11, Mmax = 5e15, virial = false, t10MF = false, xmax = 1)

This function returns a value of the 1-halo term power spectrum (in units of Mpc^3/h^3) at a specified value of the comoving wavenumber (in units of h/Mpc).

*Reference*: Equations (5,10) of Seljak, MNRAS, 318, 203 (2000)

# Arguments
- `pk`(k): a function which returns a linear matter power spectrum with the argument k being the comoving wavenumber.
   - `pk` is in units of (Mpc/h)^3.
- `z::Real`: redshift.
- `Ωm::Real`: present-day total matter density parameter.
- `k::Real`: comoving wavenumber.
   - `k` is in units of h/Mpc.

# Optional Arguments
- `Ωcb::Real=Ωm`: present-day baryon + cold dark matter density parameter.

# Optional keyword arguments
- `Mmin::Real=1e11`: minimum mass for integration, ``∫_{Mmin}^{Mmax} dM dn/dM (M/ρ)^2 |u(kR,c)|^2``.
- `Mmax::Real=5e15`: maximum mass for integration, ``∫_{Mmin}^{Mmax} dM dn/dM (M/ρ)^2 |u(kR,c)|^2``.
- `virial::Bool=false`: if `true`, use the virial overdensity `Δvir`. If `false` (the default), use `Δm=200`.
- `t10MF::Bool=false`: if `true`, use `tinker10MF` for the halo multiplicity function. If `false` (the default), use `tinker08MF`.
- `xmax::Real=1`: Fourier transform of the NFW density profile is integrated up to ``r = xmax * c * rs = xmax * r200m (if virial = false) or rvir (if virial = true)``.
"""
function onehalo(
   pk,
   z::Real,
   Ωm::Real,
   k::Real,
   Ωcb = Ωm;
   Mmin = 1e11,
   Mmax = 5e15,
   virial = false,
   t10MF = false,
   xmax = 1, 
)
   ρc = 2.775e11 # in units of h^2 M⊙/Mpc^3
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
   # %% Compute 1-halo term power spectrum, P(k) = ∫dM (dn/dM) (M/ρ)^2 |u(kR,c)|^2, in units of h^-3 Mpc3
   # c(M) relation is due to Equation (2) and "NFW, F, 0-2" of Table 1 in
   # Duffy et al., MNRAS, 390, L64 (2008)
   function dPdlnMh(lnMh)
      if virial # For Δ = Δvir(z)
         RΔh = cbrt(exp(lnMh) * 3 / 4π / (ρc * E2 * Δvir)) # in h^-1 Mpc
         A0, B0, C0 = 7.85, -0.081, -0.71 # Avir, Bvir, Cvir in Table 1 of Duffy et al.
      else # For Δ = Δm = 200
         RΔh = cbrt(exp(lnMh) * 3 / 4π / (ρc * Ωcb * Δm)) / (1 + z) # in h^-1 Mpc
         A0, B0, C0 = 10.14, -0.081, -1.01 # Amean, Bmean, Cmean in Table 1 of Duffy et al.
      end
      cnfw = A0 * (exp(lnMh) / 2e12)^B0 * (1 + z)^C0
      f = log(1 + cnfw) - cnfw / (1 + cnfw)
      c = xmax * cnfw
      # x = k*rs*(1+z), where rs is the physical NFW scale radius
      x = exp(lnk[ik]) * RΔh * (1 + z) / cnfw # IMPORTANT (1+z)
      u2 =
         (
            sin(x) * (sinint((c + 1) * x) - sinint(x)) +
            cos(x) * (cosint((c + 1) * x) - cosint(x)) -
            sin(x * c) / (c + 1) / x
         )^2 
      dPdlnMh = spl(lnMh) * exp(2 * lnMh) * u2 / f^2
   end
   res, err = hquadrature(dPdlnMh, log(Mmin), log(Mmax))
   pk_1halo = res / (ρc * Ωcb)^2 # in h^-3 Mpc^3
end
