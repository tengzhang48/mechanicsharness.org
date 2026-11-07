"""Core machinery for the grain-boundary stress-relaxation problem.

Dimensionless evolution: dSigma/dtau = L[Sigma''']  with
  L[g](Z) = PV int_0^1 s(Z,Z') g(Z') dZ',
  s(Z,Z') = 1/(Z-Z') - 1/(Z+Z') - 2 Z' (Z-Z')/(Z+Z')^3.

Eigenmodes: lambda f = -L[f'''],  f(0)=0, f'(1)=f''(1)=0.

Quadrature: substitute Z' = 1 - y^2  (y in [0,1]); then
  L[g](Z) = int_0^1 s(Z,1-y^2) G(y) dy,  G(y) = 2 y g(1-y^2).
Split s = 1/(y^2-yZ^2) + s_reg(Z,1-y^2), yZ = sqrt(1-Z), and subtract:
  L[g](Z) = int_0^1 [G(y)-G(yZ)]/(y^2-yZ^2) dy
            + G(yZ) * I0(yZ)
            + int_0^1 s_reg(Z,1-y^2) G(y) dy,
  I0(yZ) = int_0^1 dy/(y^2-yZ^2) = ln((1-yZ)/(1+yZ))/(2 yZ)  (-> -1 as yZ->0),
  s_reg(Z,Z') = -1/(Z+Z') - 2 Z'(Z-Z')/(Z+Z')^3.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import chebyshev as C


def s_reg(Z, Zp):
    return -1.0/(Z+Zp) - 2.0*Zp*(Z-Zp)/(Z+Zp)**3


def I0(yZ):
    yZ = np.asarray(yZ, dtype=float)
    out = np.empty_like(yZ)
    small = np.abs(yZ) < 1e-6
    ys = yZ[small]
    # series: ln((1-y)/(1+y))/(2y) = -(1 + y^2/3 + y^4/5 + ...)
    out[small] = -(1.0 + ys**2/3.0 + ys**4/5.0 + ys**6/7.0)
    yb = yZ[~small]
    out[~small] = np.log((1.0-yb)/(1.0+yb))/(2.0*yb)
    return out


class LOperator:
    """Applies L to functions given through G(y)=2y g(1-y^2)."""
    def __init__(self, M=512):
        yg, wg = leggauss(M)
        self.y = 0.5*(yg+1.0)      # quadrature nodes in [0,1]
        self.w = 0.5*wg
        self.M = M

    def apply(self, Z, Gfun, Gpfun):
        """L[g](Z) for Z scalar or array; Gfun(y)->G, Gpfun(y)->G'(y)."""
        Z = np.atleast_1d(np.asarray(Z, dtype=float))
        y = self.y; w = self.w
        Zp = 1.0 - y**2                       # Z' at nodes
        Gy = Gfun(y)                          # (npts,)
        res = np.empty_like(Z)
        for k, Zk in enumerate(Z):
            yZ = np.sqrt(max(1.0-Zk, 0.0))
            GZ = Gfun(np.array([yZ]))[0]
            den = (y-yZ)*(y+yZ)
            diff = (Gy - GZ)
            # patch removable singularity
            bad = np.abs(den) < 1e-13
            val = diff/np.where(bad, 1.0, den)
            if np.any(bad):
                GpZ = Gpfun(np.array([yZ]))[0]
                if yZ > 1e-14:
                    val[bad] = GpZ/(2.0*yZ)
                else:
                    # yZ=0: limit of [G(y)-G(0)]/y^2 as y->0 = G''(0)/2; use local quadratic
                    # our G are even-ish; approximate with Gpfun unavailable -> use G'(small)
                    val[bad] = GpZ/(2.0*max(yZ,1e-300)) if yZ>0 else 0.0
                    # better: [G(y)-G(0)]/y^2 at small y ~ (G(y)-G(0))/y^2 computed directly; no bad points since y>0
            part1 = np.sum(w*val)
            part2 = GZ*I0(np.array([yZ]))[0]
            part3 = np.sum(w*s_reg(Zk, Zp)*Gy)
            res[k] = part1 + part2 + part3
        return res


def S0_analytic(Z):
    """PV int_0^1 s(Z,Z') dZ'  (closed form)."""
    Z = np.asarray(Z, dtype=float)
    return (np.log((1.0+Z)/(1.0-Z)) + 6.0*Z/(1.0+Z)
            - 2.0*Z**2/(1.0+Z)**2 - 4.0)


def chebT(j, x, m=0):
    """m-th derivative of T_j at x (array)."""
    c = np.zeros(j+1); c[j] = 1.0
    for _ in range(m):
        c = C.chebder(c)
    return C.chebval(x, c)
