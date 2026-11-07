"""Reproduce the numbers marked "ours" in the keys of problems A, B, D and E (log P344).
Run: python3 reproduce_keys.py"""
import numpy as np, scipy.linalg as sl
from math import pi, sqrt, tan
from scipy.integrate import quad
# Excerpt for the grain-boundary wedge case, made on 7 November 2026: kept, unchanged, the reproduction of problem E's numbers; left out: those of problems A, B and D.

print("E (Gao et al. 1999): single-wedge eigenproblem, the paper's Gauss-Chebyshev collocation")
S = lambda x, xi: 1 / (x - xi) - 1 / (x + xi) - 2 * xi * (x - xi) / (x + xi) ** 3
Kk = lambda x, xi: np.where(x <= xi, 0.5 * xi ** 2 - 0.5 * (x - xi) ** 2, 0.5 * xi ** 2)
for N in (60, 120):
    n = 2 * N + 1
    xi = np.cos((2 * np.arange(1, N + 1) - 1) * pi / (2 * n)); x = np.cos(np.arange(1, N + 1) * pi / n)
    X, XI = np.meshgrid(x, xi, indexing="ij")
    w, Vv = sl.eig(-S(X, XI), Kk(X, XI))
    ok = np.isfinite(w) & (w.real > 0); order = np.argsort(np.where(ok, w.real, np.inf))
    R = Vv[:, order[0]].real; f = lambda xx: np.sum(pi / n * Kk(xx, xi) * R)
    Rk = np.linalg.solve(S(X, XI) * pi / n / (2 * pi), np.ones(N))
    c = np.polynomial.chebyshev.chebfit(np.concatenate([xi, -xi]), np.concatenate([Rk, Rk]), 40)
    print(f"  N={N}: lambda1 {w[order[0]].real:.6f} lambda2 {w[order[1]].real:.4f} "
          f"f1(0.5)/f1(1) {f(0.5) / f(1.0):.4f} K_inf/(s0 sqrt(pi h)) {abs(np.polynomial.chebyshev.chebval(1.0, c)) / 2:.4f}")
fB = lambda x: -0.25 * (6 * np.sqrt(1 - x * x) - 2 * (1 - x * x) ** 1.5 + 6 * x * np.arcsin(x) - 3 * pi * x - 4)
print("  Eq. B8: 3 pi / int f1 =", 3 * pi / quad(fB, 0, 1)[0])
