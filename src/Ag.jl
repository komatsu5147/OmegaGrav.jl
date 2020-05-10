function Ag(c)
   y = c - 5
   Ag = 1.0202 + 5.376e-2 * y - 1.842e-3 * y^2 + 9.489e-5 * y^3 - 2.352e-6 * y^4
end

function poly3Ag(c)
   y = c - 5
   Ag = 1.0184 + 5.389e-2 * y - 1.514e-3 * y^2 + 3.845e-5 * y^3
end

function exactAg(c)
   f = log(1 + c) - c / (1 + c)
   u2(x) =
      (
         sin(x / c) * (sinint((1 / c + 1) * x) - sinint(x / c)) +
         cos(x / c) * (cosint((1 / c + 1) * x) - cosint(x / c)) -
         sin(x) / (1 / c + 1) / x
      )^2
   res, err = hquadrature(u2, 0, 100)
   Ag = res / π / f^2
end
