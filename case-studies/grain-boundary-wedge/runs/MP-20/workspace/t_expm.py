"""Late-time decay from the *discrete* semigroup exp(A tau), which is robust:
no stiff ODE integration needed, we just apply the matrix exponential to the
initial vector and read the asymptotic slope.

Only the physically relevant (slowly decaying) part matters; the huge negative
eigenvalues contribute 3.6e-7 instantly and never affect the late-time slope.
"""

import numpy as np, time
from scipy.linalg import expm
from numpy.polynomial.legendre import leggauss
import gb3

xq, wq = leggauss(1500); xq = 0.5 * (xq + 1); wq = 0.5 * wq


def d2_dir(n):
    h = 1.0 / n; D2 = np.zeros((n + 1, n + 1))
    for i in range(1, n):
        D2[i, i] = -2.0; D2[i, i - 1] = 1.0; D2[i, i + 1] = 1.0
    D2[n, n] = -2.0; D2[n, n - 1] = 2.0
    D2[0, 0] = 1.0
    return D2 / h ** 2


for n in (12, 16, 20, 24, 28, 32):
    sn = np.arange(n + 1) * (1.0 / n)
    E = gb3.elasticity_matrix(sn, xq, wq)
    A = -(E @ d2_dir(n))
    idx = np.arange(1, n + 1)
    Ai = A[np.ix_(idx, idx)]
    evr, V = np.linalg.eig(Ai)
    real_neg = evr.real[(np.abs(evr.imag) < 1e-2 * np.abs(evr) + 1e-12) & (evr.real < 0)]
    pos = -np.sort(real_neg)[::-1]
    pos = pos[pos > 1e-3]
    ev = np.where(evr.real < 0, evr, 0.0)          # kill non-decaying modes
    g0 = np.ones(n)
    c = np.linalg.solve(V, g0)
    k = np.argmin(np.abs(sn - 0.5)) - 1
    Vk = V[k]
    taus = np.array([1.0, 1.5, 2.0, 2.5, 3.0])
    vals = np.array([np.sum(c * Vk * np.exp(ev * t)) for t in taus]).real
    m = np.abs(vals) > 1e-11
    slope = np.polyfit(taus[m], np.log(np.abs(vals[m])), 1)[0]
    print('n=%2d  lam1=%9.5f lam2=%9.4f  slope(tau)=%9.5f' % (
        n, pos[0], pos[1] if len(pos) > 1 else np.nan, -slope))