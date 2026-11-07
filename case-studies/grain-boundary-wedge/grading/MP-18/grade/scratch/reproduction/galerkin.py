import numpy as np
from numpy.polynomial.legendre import leggauss

# =============== Tailored basis Galerkin for  L s = -(1/pi) PV Int S s''' ======
# We work with the PHYSICAL kernel S_phys (decaying branch), already includes
# the 1/pi when we call M_phys.  The continuum operator on s is
#       (Ls)(x) = -(1/pi) PV Int_0^1 S_phys(x,y) s'''(y) dy.
# We seek eigenpairs L s = -lambda s (decay rate lambda>0).
#
# Basis: functions satisfying s(0)=0, s'(1)=s''(1)=0, including a x ln x term.
#   phi_0(x) = x (1 - x/2)   : s(0)=0, s'(1)=0, s''=-1 (not 0)... fix later
# We instead use simple, robust basis and enforce BCs by construction:
#   g1(x) = x^2 (3 - 2x) - 1 ... too ad hoc.
#
# Practical approach: collocation with the *integral* form solved directly.
# For a trial s with BCs, compute s''' analytically, then apply the PV integral
# with high-order quadrature to get Ls, and solve the generalized eigenproblem
# by least squares on many collocation points.

def basis_funcs(x):
    """Return list of basis functions and their 3rd derivatives.
    Basis chosen to satisfy s(0)=0, s'(1)=s''(1)=0 exactly and to include the
    x*ln x-like behavior through the x^2 ln x function."""
    x = np.asarray(x, float)
    # Use the substitution t = x ; build later smooth functions times x to vanish at 0.
    # We'll construct BC-satisfying functions:
    funcs = []
    d3s = []
    # f1 = x * ln x  (->0 at 0), but need BCs at 1. We'll combine.
    # Base set: {(x - 1)^2 * p(x)} gives s'(1)=s''(1)=0? (x-1)^2: s'=2(x-1)->0, s''=2 !=0.
    # Need s'(1)=s''(1)=0 -> factor (x-1)^3. Also s(0)=0 -> factor x (or x^k).
    # So use phi_k = x * (x-1)^3 * q_k(x), q_k polynomials. Plus x*ln x type.
    for k in range(0, 8):
        # phi_k = x*(x-1)^3 * x^k
        funcs.append(lambda xx, k=k: xx * (xx - 1) ** 3 * xx ** k)
    # log-type functions satisfying BCs: x*(x-1)^3*ln x  (ln x singular but integrable)
    for k in range(0, 6):
        funcs.append(lambda xx, k=k: xx * (xx - 1) ** 3 * xx ** k * np.log(xx + 1e-300))
    return funcs


def d3_of(f, x, h=1e-4):
    """Numerical 3rd derivative of f at x (central)."""
    return (f(x + 2 * h) - 2 * f(x + h) + 2 * f(x - h) - f(x - 2 * h)) / (2 * h ** 3)


if __name__ == "__main__":
    print("basis defined (validation in solver)")