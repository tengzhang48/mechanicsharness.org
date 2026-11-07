"""Time evolution in the OPENING variable w, which is smooth.

    g = 1 - E w ,    w_tau = D2 g = D2(1 - E w) = -D2 E w  (+0 since D2 1 = 0
                                                            only if g=1 satisfies BCs; it does NOT)
Careful: D2 1 != 0 with the Dirichlet/Neumann BCs, because g(0)=0 is imposed.
So the correct semi-discrete system is

    w_tau = D2 (1 - E w) .

We integrate this with an implicit (backward Euler) scheme, which is stable,
and monitor the decay of g(s=0.5, tau).  The asymptotic slope is lambda_1.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss
import gb3

xq, wq = leggauss(1500); xq = 0.5 * (xq + 1); wq = 0.5 * wq


def d2_dir(s):
    n = len(s) - 1
    D2 = np.zeros((n + 1, n + 1))
    for i in range(1, n):
        hm = s[i] - s[i - 1]; hp = s[i + 1] - s[i]
        D2[i, i - 1] = 2.0 / (hm * (hm + hp))
        D2[i, i + 1] = 2.0 / (hp * (hm + hp))
        D2[i, i] = -(D2[i, i - 1] + D2[i, i + 1])
    hm = s[n] - s[n - 1]
    D2[n, n] = 2.0 / hm ** 2
    D2[n, n - 1] = -2.0 / hm ** 2
    D2[0, 0] = 1.0
    return D2


n = 40
s = np.arange(n + 1) * (1.0 / n)
E = gb3.elasticity_matrix(s, xq, wq)
D2 = d2_dir(s)
# g = 1 - E w  with g_0 = 0 pinned: over-write the s=0 row of the elasticity map
# g(s) = 1 - (E w)(s) for interior, and g(0)=0.
E[0, :] = 0.0

w = np.zeros(n + 1)
dt = 0.02
k = np.argmin(np.abs(s - 0.5))
hist = []
# implicit: w_new = w + dt*D2*(1 - E w_new)
M = np.eye(n + 1) + dt * (D2 @ E)
for step in range(1, 50001):
    rhs = w + dt * (D2 @ np.ones(n + 1))
    w = np.linalg.solve(M, rhs)
    if step % 2500 == 0:
        g = 1.0 - E @ w
        g[0] = 0.0
        hist.append((step * dt, g[k]))
        print('tau=%8.3f  g(0.5)=%+.6e' % (step * dt, g[k]))

tt = np.array([h[0] for h in hist]); gg = np.abs(np.array([h[1] for h in hist]))
m = (tt > 300.0) & (gg > 1e-13)
if np.sum(m) > 3:
    p = np.polyfit(tt[m], np.log(gg[m]), 1)
    print('late-time slope -> lambda_1 = %.5f' % (-p[0]))