"""
Corrected shooting: only ONE axis BC (C(0)=0 => phi'(0)=0).  Free ICs:
phi(0)=1 (scale fix), phi'(0)=0, and (phi''(0), phi'''(0)) = (u1,u2) unknown.
Solve 2x2: A(pi)=0, C(pi)=0.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import hrr2 as H


def endvals(n, u1, u2):
    H.n = n
    H.p = 1.0/(n+1.0)
    H.k = 2.0 - H.p
    Y0 = [1.0, 0.0, u1, u2]
    s = solve_ivp(H.rhs, [0.0, np.pi], Y0, rtol=1e-10, atol=1e-12,
                  method='Radau', max_step=0.01)
    if not s.success:
        return np.array([1e3, 1e3])
    phi, d1, d2, d3 = s.y[:, -1]
    return np.array([H.k*phi + d2, -(H.k-1.0)*d1])


def solve(n, guess):
    H.n = n
    H.p = 1.0/(n+1.0)
    H.k = 2.0 - H.p
    sol, info, ier, msg = fsolve(lambda x: endvals(n, x[0], x[1]), guess,
                                 xtol=1e-11, full_output=True)
    return sol, ier, msg


for n, g in [(1.5, [1.0, 0.0]), (3.0, [0.2, 0.0]), (6.5, [0.0, 0.0])]:
    s, ier, msg = solve(n, g)
    r = endvals(n, s[0], s[1])
    print("n=%.2f  phi''(0)=%.6f phi'''(0)=%.6f  A(pi)=%.2e C(pi)=%.2e ier=%d"
          % (n, s[0], s[1], r[0], r[1], ier))