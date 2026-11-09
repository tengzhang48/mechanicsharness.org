"""
Robust shooting for the HRR Airy ODE with continuation in n.

State Y=[phi,phi',phi'',phi'''] on [0,pi].
Symmetry: phi'(0)=0, phi'''(0)=0.
Free surface: A(pi)=k phi+phi''=0  AND  C(pi)=-(k-1)phi'=0.

With 2 free ICs (phi(0), phi''(0)) and 2 free-surface conditions, the system is
homogeneous in scale -> solve for the ratio r = phi''(0)/phi(0) with phi(0)=1
using the TWO conditions. Over-determination is resolved only if the two
conditions are dependent. Numerically we will find the r that zeros A(pi) and
report C(pi) there; if C(pi)!=0 the leading self-similar form is inconsistent
and we report that finding.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import hrr2 as H


def endvals(n, t):
    H.n = n
    H.p = 1.0/(n+1.0)
    H.k = 2.0 - H.p
    Y0 = [1.0, 0.0, t, 0.0]
    s = solve_ivp(H.rhs, [0.0, np.pi], Y0, rtol=1e-11, atol=1e-13,
                  method='Radau', max_step=0.005)
    if not s.success:
        return None
    phi, d1, d2, d3 = s.y[:, -1]
    A = H.k*phi + d2
    C = -(H.k-1.0)*d1
    return A, C


def find_t_A0(n, tlo=-4, thi=4, N=400):
    ts = np.linspace(tlo, thi, N)
    prev = None
    roots = []
    for t in ts:
        v = endvals(n, t)
        if v is None:
            prev = None; continue
        A = v[0]
        if prev is not None and prev[1]*A < 0:
            try:
                r = brentq(lambda x: endvals(n, x)[0], prev[0], t, xtol=1e-12)
                roots.append(r)
            except Exception:
                pass
        prev = (t, A)
    return roots


if __name__ == "__main__":
    for n in [1.1, 1.5, 2.0, 3.0, 4.0, 6.5]:
        roots = find_t_A0(n)
        print("n=%.2f" % n, "roots t* (A(pi)=0):",
              ["%.5f" % r for r in roots])
        for r in roots:
            A, C = endvals(n, r)
            print("    t*=%.5f  A(pi)=%.2e  C(pi)=%.5f" % (r, A, C))