"""Grader audit, using only the copied submission and local computation."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent / "submission"))
import numpy as np
import scipy.linalg as sla
from numpy.polynomial.chebyshev import chebval
import solve

A, F, xc = solve.build_pencil(48, 800)
lam, V = sla.eig(A, -F)
sel = (np.abs(lam.imag) < 1e-8 * np.maximum(np.abs(lam.real), 1e-3)) & (lam.real > 1e-6)
idx = np.where(sel)[0]
order = idx[np.argsort(lam[idx].real)]
lams = lam[order].real
a = np.real(V[:, order[0]])
print("Computed first twelve rates:", lams[:12])
print("Actual first five successive ratios:", lams[1:6] / lams[:5])
n = np.arange(1, 6)
print("Ratios required by lambda_n proportional to n^4.5:", ((n + 1) / n) ** 4.5)
print("lambda_n / n^3, n=6..12:", lams[5:12] / np.arange(6, 13) ** 3)

t, wt = solve.tquad(0.0, 1.0, 2000)
xi = 1 - t**2
vq = chebval(2 * xi - 1, a)
Gq = xi / t * vq
xx = np.linspace(0.0005, 0.9995, 2000)
rr = []
for x in xx:
    Gx = x / np.sqrt(1 - x) * chebval(2 * x - 1, a)
    I1 = 2 * np.sum(wt * (t * Gq - t * Gx) / (x - 1 + t**2))
    I2 = Gx * np.log(x / (1 - x))
    I3 = 2 * np.sum(wt * solve.Rfun(x, xi) * t * Gq)
    rr.append(abs(I1 + I2 + I3 + lams[0] * solve.f_of(a, x)[0]))
rr = np.array(rr)
im = np.argmax(rr)
print("Actual independent 2000-point residual grid: max", rr[im], "at x", xx[im])

Nm = 10
xg = np.linspace(0.05, 0.95, 19)
Fm = np.array([solve.f_of(np.real(V[:, j]), xg) for j in order[:Nm]])
c, *_ = np.linalg.lstsq(Fm.T, np.ones_like(xg), rcond=None)
print("Ten-mode IC reconstruction max error on fit grid:", np.max(np.abs(Fm.T @ c - 1)))
taus = np.linspace(0.5, 1.2, 8)
wrong_rate = 2 * lams[0]
wrong_vals = c[0] * solve.f_of(a, 0.5)[0] * np.exp(-wrong_rate * taus)
print("Injected twice-too-large slow rate and fitted decay:", wrong_rate,
      np.polyfit(taus, np.log(wrong_vals), 1)[0])
for nm, ng in [(5, 19), (10, 19), (10, 199), (16, 399)]:
    grid = np.linspace(0.001, 0.999, ng) if ng > 19 else xg
    modes = [np.real(V[:, j]) for j in order[:nm]]
    mat = np.array([solve.f_of(av, grid) for av in modes])
    coeff, *_ = np.linalg.lstsq(mat.T, np.ones_like(grid), rcond=None)
    amp05 = coeff[0] * solve.f_of(modes[0], 0.5)[0]
    print("IC projection nm/ng:", nm, ng, "late mid-boundary amplitude:", amp05,
          "threshold times:", np.log(amp05 / np.array([1 / np.e, 0.05])) / lams[0])

import eig_g
print("Delivered mode_shape return:", eig_g.mode_shape(lams[0], V, xc))
from gbcore import diff_mat
arbitrary_nodes = np.array([0.0, 0.2, 0.5, 1.0])
print("General-node diff_mat max derivative error for x^2:",
      np.max(np.abs(diff_mat(arbitrary_nodes) @ arbitrary_nodes**2 - 2 * arbitrary_nodes)))
from eig import eigen_s
x, eigvals, eigvecs, positive = eigen_s(N=24, verbose=False)
print("eig.eigen_s returned eigenvalue/vector dimensions:", len(eigvals), eigvecs.shape)
from crack import crack_solve
kd, ka, _ = crack_solve(N=32, verbose=False)
print("Development crack.py density and finite-distance ahead K estimates:", kd, ka)
