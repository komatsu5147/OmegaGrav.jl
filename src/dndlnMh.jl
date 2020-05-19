"""
    dndlnMh(pk, z, Δm, Ωm, lnMh; t10MF = false)

This function returns a value of the comoving halo number density per log mass interval in units of h^3/Mpc^3 at a specified value of the log(halo mass).

# Arguments
- `pk`(k): a function which returns a linear matter power spectrum with the argument k being the comoving wavenumber.
   - `pk` is in units of (Mpc/h)^3.
- `z::Real`: redshift.
- Δm: halo overdensity with respect to the mean mass density of the Universe
- `Ωm::Real`: present-day matter density parameter. When excluding the neutrino contribution, use Ωcb (baryon + cold dark matter density parameter).
- `lnMh::Real`: log(halo mass)
   - halo mass is in units of M⊙/h.

# Optional keyword arguments
- `t10MF::Bool=false`: if `true`, use `tinker10MF` for the halo multiplicity function. If `false` (the default), use `tinker08MF`.
"""
function dndlnMh(pk, z, Δm, Ωm, lnMh; t10MF = false)
   ρc = 2.775e11 # in units of h^2 M⊙/Mpc^3
   Rh = cbrt(exp(lnMh) * 3 / 4π / ρc / Ωm) # in units of Mpc/h
   σ2 = sigma2(pk, Rh)
   dlnσ2dlnRh = Rh * dsigma2dR(pk, Rh) / σ2
   lnν = 2 * log(1.6865) - log(σ2)
   if t10MF
      MF = tinker10MF(lnν, z, Δm)
   else
      MF = tinker08MF(lnν, z, Δm)
   end
   dndlnMh = -dlnσ2dlnRh * MF / 4π / Rh^3 # in units of h^3 Mpc^-3
end
