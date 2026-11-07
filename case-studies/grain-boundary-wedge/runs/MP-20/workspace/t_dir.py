import numpy as np, time
from numpy.polynomial.legendre import leggauss
import gb3

xq, wq = leggauss(1500); xq = 0.5 * (xq + 1); wq = 0.5 * wq


def d2_dir(n):
    h = 1.0 / n; D2 = np.zeros((n + 1, n + 1))
    for i in range(1, n):
        D2[i, i] = -2.0; D2[i, i - 1] = 1.0; D2[i, i + 1] = 1.0
    D2[n, n] = -2.0; D2[n, n - 1] = 2.0      # Neumann at s=1
    D2[0, 0] = 1.0                            # Dirichlet at s=0
    return D2 / h ** 2


print('  n      lam1       lam2       lam3      |imag|max')
for n in (10, 14, 18, 22, 26, 30, 34, 40, 50):
    t = time.time()
    sn = np.arange(n + 1) * (1.0 / n)
    E = gb3.elasticity_matrix(sn, xq, wq)
    A = -(E @ d2_dir(n))
    idx = np.arange(1, n + 1)
    Ai = A[np.ix_(idx, idx)]
    ev = np.linalg.eigvals(Ai)
    evr = np.sort(ev.real[np.abs(ev.imag) < 1e-2 * np.abs(ev)])
    pos = evr[evr > 1e-3]
    s3 = pos[2] if len(pos) > 2 else float('nan')
    print('%3d  %9.5f %9.4f %9.4f   %.1e  (%.1fs)' % (
        n, pos[0], pos[1], s3, np.max(np.abs(ev.imag)), time.time() - t))