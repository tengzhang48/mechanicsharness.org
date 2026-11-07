import numpy as np
from numpy.polynomial import chebyshev as C
from numpy.polynomial.legendre import leggauss

# =============== Spectral Galerkin with exact 3rd derivatives ===============
# Represent  s(x) = sum_k a_k T_k(t),  t = 2x-1 in [-1,1].
# Impose BCs s(0)=0, s'(1)=0, s''(1)=0 by choosing a basis of functions that
# satisfy them exactly (build nullspace of the 3 BC linear functionals in the
# Chebyshev coefficient space).  Then s''' is exact from the series.
#
# Operator: (Ls)(x) = -(1/pi) PV Int_0^1 S_phys(x,y) s'''(y) dy.
# Galerkin: project onto the same BC-satisfying space:
#    <psi_j, L psi_i> inner products when psi are polynomials -> but kernel is
#    singular.  Use TEST functions that are polynomials and compute the double
#    integral by swapping: Int psi_j(x) Int S_phys(x,y) s'''(y) dy dx.
#    The inner x-integral of S_phys(x,y) over x with a polynomial psi_j(x) can be
#    done with the SAME antiderivative machinery (roles swapped), giving a
#    regular function of y times s'''(y) -> then a 1D quadrature in y.


def cheb_basis_with_BCs(Nc):
    """Return coefficient vectors (in Chebyshev basis) for functions satisfying
    s(0)=0, s'(1)=0, s''(1)=0.  We restrict to the nullspace of the 3 linear
    functionals acting on Chebyshev coefficients up to degree Nc."""
    deg = np.arange(Nc + 1)

    def evals(cvec, x):
        return C.chebval(2 * x - 1, cvec)

    # Build functional matrices via exact coefficient vectors of T_k and its derivs
    # at x=0 (t=-1) and x=1 (t=1).
    # s(0): T_k(-1) = (-1)^k
    A0 = np.array([(-1.0) ** k for k in deg])
    # d/dx T_k = 2 dT_k/dt ; dT_k/dt at t=1 = k^2
    def dTdt_at(c, t):
        # numerical exactness: use C.chebder
        dc = C.chebder(c)
        return C.chebval(t, dc)
    # first derivative at x=1
    A1 = np.array([2 * dTdt_at(np.eye(Nc + 1)[k], 1.0) for k in deg])
    # second derivative at x=1
    A2 = np.array([4 * dTdt_at(C.chebder(np.eye(Nc + 1)[k]), 1.0) for k in deg])
    Mbc = np.vstack([A0, A1, A2])       # 3 x (Nc+1)
    # nullspace
    u, s, vt = np.linalg.svd(Mbc)
    # vt rows with s ~ 0 give the nullspace
    ns = vt[3:]                          # (Nc+1-3) x (Nc+1)
    return ns


if __name__ == "__main__":
    ns = cheb_basis_with_BCs(12)
    print("basis size:", ns.shape)
    # verify BCs
    x0, x1 = 0.0, 1.0
    for i in range(ns.shape[0]):
        c = ns[i]
        v0 = C.chebval(2 * x0 - 1, c)
        d1 = 2 * C.chebval(2 * x1 - 1, C.chebder(c))
        d2 = 4 * C.chebval(2 * x1 - 1, C.chebder(C.chebder(c)))
        if i < 4:
            print("basis %d: s(0)=%.2e s'(1)=%.2e s''(1)=%.2e" % (i, v0, d1, d2))