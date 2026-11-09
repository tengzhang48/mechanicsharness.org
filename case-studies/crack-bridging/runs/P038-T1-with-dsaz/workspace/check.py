"""Independent check of the critical state, Solver B.

Semi-infinite crack in a K-field (front-scale model).  Zone 0<x<X with
uniform closure sigma0; tip at x=0.  Governing equations from zone.py:
  K_tip = K0 - sqrt(2/pi) * (2 sigma0 sqrt(X))            [Dugdale-consistent]
  delta(x) = (8 K0/E') sqrt(x/2pi)
             - (4/(pi E')) int_0^X ln|(sqrt x+sqrt xi)/(sqrt x-sqrt xi)| p dxi
Here the closure is uniform, so the integral is closed form:
  J(x) = int_0^X ln|(sqrt x+sqrt xi)/(sqrt x-sqrt xi)| dxi.
Front condition: K_tip = K0 (with E'=1, sigma0=1).  Then
  sqrt(2/pi)*2*sqrt(X) = 0 ... no: K_tip=Kfront=K0 means closure only
  reduces... Actually we impose the *intrinsic* front: K_tip = K0intrinsic
  = sqrt(Gamma0) and the applied front K = K0intrinsic + closure drop.
This model NEEDS the exterior (finite crack) to fix the relation between
the applied K and the geometry; it is used here only for the *ratio*
delta_t / (K0 sqrt X) scaling as a cross-check of the zone profile shape,
not as an absolute answer.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss


def J_integral(x, X, M=400):
    """int_0^X ln|(sqrt x + sqrt xi)/(sqrt x - sqrt xi)| dxi (closed away
    from the log singularity by Gauss split)."""
    if x <= 1e-15:
        return 0.0
    t, w = leggauss(M)
    total = 0.0
    for x0, x1 in ((0.0, x), (x, X)):
        if x1 - x0 < 1e-300:
            continue
        xi = 0.5 * (x1 - x0) * t + 0.5 * (x1 + x0)
        wi = 0.5 * (x1 - x0) * w
        ker = np.log(np.abs((np.sqrt(x) + np.sqrt(xi)) /
                            (np.sqrt(x) - np.sqrt(xi))))
        ker[~np.isfinite(ker)] = 0.0
        total += np.sum(wi * ker)
    return total


def delta_profile(X, Kfront, sigma0=1.0, E=1.0, n=60):
    """delta(x) for uniform closure sigma0 on 0<x<X, applied front K=Kfront."""
    xs = np.linspace(0, X, n)
    d = np.array([(8.0 * Kfront / E) * np.sqrt(x / (2.0 * np.pi))
                  - (4.0 / (np.pi * E)) * sigma0 * J_integral(x, X, M=160)
                  for x in xs])
    return xs, d


if __name__ == '__main__':
    # Normalized: at the critical state E'=1, sigma0=1, delta_c=1.
    # Ln* ~ 0.52, an ~ 0.3025.  The front K (local) = K0n + closure drop.
    # closure drop = sqrt(2/pi)*2*sqrt(X); X = Ln (zone length per end).
    Ln = 0.520121
    K0n = np.sqrt(Gamma0 / (sigma0 * delta_c)) if False else None
    from solve import Gamma0, sigma0 as s0, delta_c as dc
    K0n = np.sqrt(Gamma0 / (s0 * dc))
    X = Ln
    Kfront = K0n + np.sqrt(2.0 / np.pi) * 2.0 * np.sqrt(X)
    xs, d = delta_profile(X, Kfront)
    # delta at the root x=X (trailing edge) vs at the front x=0
    print(f"front-model: Kfront={Kfront:.4f}, K0n={K0n:.4f}, X={X:.4f}")
    print(f"  delta(root x=X) = {d[-1]:.4f}")
    print(f"  delta(front x=0) -> 0 by construction, give delta(X*0.5)={np.interp(0.5*X,xs,d):.4f}")