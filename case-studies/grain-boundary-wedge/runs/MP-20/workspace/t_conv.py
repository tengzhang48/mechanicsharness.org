import numpy as np, gb3, time
from numpy.polynomial.legendre import leggauss

xq, wq = leggauss(1500); xq = 0.5 * (xq + 1); wq = 0.5 * wq


def d2_tri(n):
    h = 1.0 / n; D2 = np.zeros((n, n))
    for i in range(n):
        D2[i, i] = -2.0
        if i > 0: D2[i, i - 1] = 1.0
        if i < n - 1: D2[i, i + 1] = 1.0
    D2[0, 0] = -1.0; D2[n - 1, n - 1] = -1.0
    return D2 / h ** 2


print(' n     lam1       lam2       lam3        |imag|')
for n in (8, 12, 16, 20, 26, 34, 44, 60, 80, 100, 130, 170, 220):
    t = time.time()
    sn = (np.arange(n) + 0.5) / n
    E = gb3.elasticity_matrix(sn, xq, wq)
    A = -(E @ d2_tri(n))
    ev = np.linalg.eigvals(A)
    real = ev[np.abs(ev.imag) < 1e-2 * np.abs(ev)]
    evr = np.sort(real.real)
    pos = evr[evr > 1.0]
    s3 = pos[2] if len(pos) > 2 else float('nan')
    print('%3d  %10.5f %10.5f %10.5f  |imag|=%.1e  (%.1fs)' % (
        n, pos[0], pos[1], s3, np.max(np.abs(ev.imag)), time.time() - t))