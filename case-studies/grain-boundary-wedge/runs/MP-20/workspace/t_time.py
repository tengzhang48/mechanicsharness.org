"""Independent confirmation of lambda_1 by direct time integration, with the
surface boundary condition g(0)=0 enforced (that is the physics: the traction
is pinned to zero at the free surface, and that drives the relaxation).

Grid: s_i = i*h, i = 0..n (h=1/n), with g_0 = 0 fixed.  Unknowns are interior
nodal values of g.  Diffusion operator D2 is the three-point second difference
with Dirichlet at 0 and Neumann (g_s(1)=0) at 1.  System:  g_tau = -E D2 g.
"""

import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss
import gb3

xq, wq = leggauss(1500); xq = 0.5 * (xq + 1); wq = 0.5 * wq


def build(n):
    sn = np.arange(n + 1) * (1.0 / n)
    E = gb3.elasticity_matrix(sn, xq, wq)
    h = 1.0 / n
    D2 = np.zeros((n + 1, n + 1))
    for i in range(1, n):
        D2[i, i] = -2.0; D2[i, i - 1] = 1.0; D2[i, i + 1] = 1.0
    D2[n, n] = -2.0; D2[n, n - 1] = 2.0
    D2[0, 0] = 1.0
    D2 /= h * h
    return sn, E, D2


n = 30
sn, E, D2 = build(n)
A = -(E @ D2)
idx = np.arange(1, n + 1)
Ai = A[np.ix_(idx, idx)]
ev = np.linalg.eigvals(Ai)
evr = np.sort(ev.real[np.abs(ev.imag) < 1e-3 * np.abs(ev)])
print('smallest real eigenvalues (n=%d):' % n, evr[:6])

rhs = lambda tau, gi: Ai @ gi
g0 = np.ones(n)
sol = solve_ivp(rhs, [0, 5.0], g0, method='BDF', rtol=1e-10, atol=1e-14,
                dense_output=True, max_step=0.05)
k = np.argmin(np.abs(sn - 0.5)) - 1
print('tau, g(0.5):')
for t in np.linspace(0, 5.0, 26):
    print('  %7.3f  %+.6e' % (t, sol.sol(t)[k]))
tt = np.linspace(2.5, 5.0, 200)
gg = np.abs(np.array([sol.sol(t)[k] for t in tt]))
good = gg > 1e-14
p = np.polyfit(tt[good], np.log(gg[good]), 1)
print('fitted lambda_1 = %.4f' % (-p[0]))