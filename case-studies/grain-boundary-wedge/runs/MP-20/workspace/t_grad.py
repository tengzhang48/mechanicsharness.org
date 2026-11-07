"""
Graded-grid discretisation of the operator A = -E D2 with the surface
condition g(0)=0, designed to resolve the logarithmic layer at s = 0.

Grid: s_i = (i/n)^q with q >= 2, i = 0..n (clustered at the surface).
Diffusion uses the standard three-point second difference on the non-uniform
grid with Dirichlet at s=0 and Neumann at s=1.  The elasticity matrix is
evaluated with the validated p.v. quadrature.  We monitor the two smallest
positive (decay) rates as q and n are refined.
"""

import numpy as np, time
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
    D2[n, n] = 2.0 / hm ** 2          # Neumann: reflect
    D2[n, n - 1] = -2.0 / hm ** 2
    D2[0, 0] = 1.0
    return D2


print('  n    q      lam1        lam2')
for q in (1.0, 1.5, 2.0, 2.5, 3.0):
    for n in (16, 24, 32, 48):
        t = time.time()
        s = (np.arange(n + 1) / n) ** q
        E = gb3.elasticity_matrix(s, xq, wq)
        A = -(E @ d2_graded(s))
        idx = np.arange(1, n + 1)
        Ai = A[np.ix_(idx, idx)]
        ev = np.linalg.eigvals(Ai)
        real = ev[np.abs(ev.imag) < 1e-2 * np.abs(ev)].real
        pos = np.sort(real[real > 1e-3])
        print('%3d  %.1f  %10.6f  %10.4f   (%.1fs)' % (
            n, q, pos[0], pos[1] if len(pos) > 1 else np.nan, time.time() - t))