import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import hrr2 as H


def endvals(n, t):
    H.n = n
    H.p = 1.0/(n+1.0)
    H.k = 2.0 - H.p
    Y0 = [1.0, 0.0, t, 0.0]
    s = solve_ivp(H.rhs, [0.0, np.pi], Y0, rtol=1e-9, atol=1e-11,
                  method='Radau', max_step=0.02)
    if not s.success:
        return None
    phi, d1, d2, d3 = s.y[:, -1]
    return H.k*phi + d2, -(H.k-1.0)*d1   # A(pi), C(pi)


for n in [1.5, 3.0, 6.5]:
    ts = np.linspace(-3, 3, 60)
    vals = [(t, endvals(n, t)) for t in ts]
    roots = []
    prev = None
    for t, v in vals:
        if v is None:
            prev = None; continue
        A = v[0]
        if prev is not None and prev[1]*A < 0:
            try:
                r = brentq(lambda x: endvals(n, x)[0], prev[0], t, xtol=1e-10)
                roots.append(r)
            except Exception as e:
                pass
        prev = (t, A)
    print("n=%.2f  roots(A(pi)=0):" % n, ["%.5f" % r for r in roots])
    for r in roots:
        A, C = endvals(n, r)
        print("     t*=%.5f A=%.2e C(pi)=%.5f" % (r, A, C))