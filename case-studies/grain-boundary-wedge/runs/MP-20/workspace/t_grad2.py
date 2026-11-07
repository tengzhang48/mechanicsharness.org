"""lambda_1 on strongly graded grids, moderate n, looking for convergence."""
import numpy as np
from numpy.polynomial.legendre import leggauss
import gb3

xq, wq = leggauss(2000); xq = 0.5 * (xq + 1); wq = 0.5 * wq


def d2_graded(s):
    n = len(s) - 1
    D2 = np.zeros((n + 1, n + 1))
    for i in range(1, n):
        hm = s[i] - s[i - 1]; hp = s[i + 1] - s[i]
        D2[i, i - 1] = 2.0 / (hm * (hm + hp))
        D2[i, i + 1] = 2.0 / (hp * (hm + hp))
        D2[i, i] = -(D2[i, i - 1] + D2[i, i + 1])
    hm = s[n] - s[n - 1]
    D2[n, n] = 2.0 / hm ** 2; D2[n, n - 1] = -2.0 / hm ** 2
    D2[0, 0] = 1.0
    return D2


print(' n   q       lam1        lam2      lam3')
for q in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
    for n in (12, 16, 20, 24):
        s = (np.arange(n + 1) / n) ** q
        E = gb3.elasticity_matrix(s, xq, wq)
        A = -(E @ d2_graded(s))
        idx = np.arange(1, n + 1)
        Ai = A[np.ix_(idx, idx)]
        ev = np.linalg.eigvals(Ai)
        real = ev[np.abs(ev.imag) < 1e-3 * np.abs(ev)].real
        pos = np.sort(real[real > 1e-4])
        print('%3d  %.1f %11.6f %11.4f %11.4f' % (
            n, q, pos[0], pos[1] if len(pos) > 1 else np.nan,
            pos[2] if len(pos) > 2 else np.nan))