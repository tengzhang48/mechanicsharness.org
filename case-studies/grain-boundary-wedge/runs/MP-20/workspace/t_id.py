"""Identify the physical slowest mode on uniform grids and extrapolate.

The suspect mode: the one whose eigenvector is smooth and O(1) in the interior
(not a surface spike).  We track it over n and Richardson-extrapolate.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
import gb3

xq, wq = leggauss(2000); xq = 0.5 * (xq + 1); wq = 0.5 * wq


def d2_uniform(n):
    h = 1.0 / n; D2 = np.zeros((n + 1, n + 1))
    for i in range(1, n):
        D2[i, i] = -2.0; D2[i, i - 1] = 1.0; D2[i, i + 1] = 1.0
    D2[n, n] = -2.0; D2[n, n - 1] = 2.0
    D2[0, 0] = 1.0
    return D2 / h ** 2


res = []
for n in (8, 10, 12, 14, 16, 18, 20, 22, 24):
    s = np.arange(n + 1) * (1.0 / n)
    E = gb3.elasticity_matrix(s, xq, wq)
    A = -(E @ d2_uniform(n))
    idx = np.arange(1, n + 1)
    Ai = A[np.ix_(idx, idx)]
    ev, V = np.linalg.eig(Ai)
    real = np.abs(ev.imag) < 1e-3 * np.abs(ev)
    # physical mode: smoothest eigenvector with small |lambda|
    best = None
    for j in np.where(real & (ev.real > 1e-4))[0]:
        v = V[:, j].real; v = v / np.max(np.abs(v))
        rough = np.sum(np.abs(np.diff(v)))   # roughness
        if best is None or (ev.real[j] < 3.0 and rough < 2.0):
            if best is None or ev.real[j] < best[0]:
                best = (ev.real[j], rough)
    inner = [ev.real[j] for j in np.where(real & (ev.real > 1e-4))[0]]
    inner.sort()
    print('n=%2d  first few positive: %s' % (n, np.array2string(np.array(inner[:4]), precision=5)))
    res.append((n, inner[0] if inner else np.nan))

res = np.array(res)
print('\nRichardson on n (lam1 vs 1/n):')
p = np.polyfit(1.0 / res[:, 0], res[:, 1], 1)
print('  fit lam1 = %.5f + %.5f*(1/n)  -> intercept %.5f' % (p[1], p[0], p[1]))