"""Problem E key values by the paper's Gauss-Chebyshev collocation (Gao et al. 1999, App. A, Eqs. A5-A8),
extended to N = 240, 480 (audit 2026-10-25). Same scheme as reproduce_keys.py and the bundle's
grain-boundary_check.py; adds larger N and Richardson extrapolation (error ~ 1/N assumed, checked by the
successive differences). Also the mode-1 amplitude c1 (not requested by the problem) and the half-time
t_1/2 at z = h, as a candidate well-defined replacement for part (c).
Usage: <volume>/venvs/coupfe-eda/bin/python -I key_recompute.py > key_recompute.out"""
import math
import numpy as np, scipy.linalg as sl
pi = math.pi
S = lambda x, xi: 1 / (x - xi) - 1 / (x + xi) - 2 * xi * (x - xi) / (x + xi) ** 3
Kk = lambda x, xi: np.where(x <= xi, 0.5 * xi ** 2 - 0.5 * (x - xi) ** 2, 0.5 * xi ** 2)
rows = {}
for N in (30, 60, 120, 240, 480):
    n = 2 * N + 1
    th = (2 * np.arange(1, N + 1) - 1) * pi / (2 * n); xi = np.cos(th); x = np.cos(np.arange(1, N + 1) * pi / n)
    X, XI = np.meshgrid(x, xi, indexing="ij")
    A, Bm = -(pi / n) * S(X, XI), (pi / n) * Kk(X, XI)
    w, V = sl.eig(A, Bm)
    ok = np.isfinite(w) & (np.abs(w.imag) < 1e-8 * np.abs(w.real)) & (w.real > 0)
    order = np.argsort(np.where(ok, w.real, np.inf))
    v1 = V[:, order[0]]
    f = lambda xx, v: np.sum(pi / n * Kk(xx, xi) * v)
    shape = (f(0.5, v1) / f(1.0, v1)).real
    # traction-free limit: (pi/n) S R = 1 ; K/(sigma0 sqrt(pi h)) = pi |R(1)| via even-Chebyshev interpolation of R
    Rk = np.linalg.solve((pi / n) * S(X, XI), np.ones(N))
    T = np.cos(2 * np.outer(th, np.arange(N))); a = np.linalg.solve(T, Rk)
    Kinf = pi * abs(np.sum(a))
    # c1: expand the discrete uniform initial traction on the discrete modes (f(1) = 1 normalization)
    r0 = np.linalg.solve(Bm, np.ones(N)); c = np.linalg.solve(V, r0.astype(complex))
    c1 = (c[order[0]] * np.sum((pi / n) * 0.5 * xi ** 2 * v1)).real
    rows[N] = (w[order[0]].real, w[order[1]].real, shape, Kinf, c1)
    print(f"N={N:4d}: lambda1 {rows[N][0]:.7f} lambda2 {rows[N][1]:.5f} f1(0.5)/f1(1) {shape:.6f} K_inf {Kinf:.5f} c1 {c1:.5f}")
print("Richardson (2*X(2N) - X(N)), assuming error ~ 1/N:")
names = ["lambda1", "lambda2", "shape", "K_inf", "c1"]
for N in (120, 240):
    ext = [2 * rows[2 * N][i] - rows[N][i] for i in range(5)]
    print(f"  from N={N},{2*N}: " + "  ".join(f"{nm} {v:.6f}" for nm, v in zip(names, ext)))
l1 = 11.742668; c1 = 2 * rows[480][4] - rows[240][4]
print(f"e-folding time 1/lambda1 = {1 / l1:.5f} t0")
print(f"half-time at z = h (higher modes decay by exp(-lambda2 t) ~ e^-25 there): ln(2 c1)/lambda1 = {math.log(2 * c1) / l1:.5f} t0 (c1 = {c1:.5f})")
print(f"95% relaxation at z = h: ln(20 c1)/lambda1 = {math.log(20 * c1) / l1:.4f} t0;  ln(20)/lambda1 = {math.log(20) / l1:.4f} t0")

# ---- appended: half-time at z = h from the full discrete modal sum (all modes), N = 120, 240, 480
from scipy.optimize import brentq
print("half-time t_1/2 at z = h from the full modal sum sigma(1,t) = sum c_n f_n(1) exp(-lambda_n t):")
for N in (120, 240, 480):
    n = 2 * N + 1
    th = (2 * np.arange(1, N + 1) - 1) * pi / (2 * n); xi = np.cos(th); x = np.cos(np.arange(1, N + 1) * pi / n)
    X, XI = np.meshgrid(x, xi, indexing="ij")
    A, Bm = -(pi / n) * S(X, XI), (pi / n) * Kk(X, XI)
    w, V = sl.eig(A, Bm)
    good = np.isfinite(w) & (w.real > 0)
    r0 = np.linalg.solve(Bm, np.ones(N)); c = np.linalg.solve(V, r0.astype(complex))
    f1n = (pi / n) * 0.5 * xi ** 2 @ V                     # f_n(1) for every mode
    sig1 = lambda t: float(np.real(np.sum((c * f1n * np.exp(-w * t))[good])))
    print(f"  N={N}: sigma(1, 0.001) = {sig1(0.001):.4f}, sigma(1, 0.05) = {sig1(0.05):.4f}, t_1/2 = {brentq(lambda t: sig1(t) - 0.5, 0.01, 0.5):.5f} t0")
