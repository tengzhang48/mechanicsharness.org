"""Grain-boundary diffusion wedge (multiphysics round, problem E): recompute the reference numbers.

Reference: <study-folder>/multiphysics_round/E/KEY.md, from H. Gao, L. Zhang, W. D. Nix,
C. V. Thompson and E. Arzt, Acta Mater. 47 (1999) 2865-2878 (Eqs. 13-18, Appendix A Eqs. A3-A8, Appendix B).

1. The eigenproblem  lambda f(x) = -int_0^1 S(x,xi) f'''(xi) dxi,  f(0) = f'(1) = f''(1) = 0,  written as the
   paper's Cauchy-type equation (Eq. 18) for f''' = R(xi)/sqrt(1-xi^2) and solved by its Gauss-Chebyshev
   collocation (Eqs. A5-A8), N = 30, 60, 120. R is interpolated in even Chebyshev polynomials to evaluate
   f1(x) and the opening shape -f1''(x) off the nodes.
2. The traction-free limit: int_0^1 S v dxi = 1 with v = (E*/(2 pi sigma0 h)) du/dxi = R/sqrt(1-xi^2); then
   u -> (4K/E*) sqrt(r/(2 pi)) gives K/(sigma0 sqrt(pi h)) = pi |R(1)|.
3. Appendix B: f1'' = -3 sqrt(1-x^2) (Eq. B6), lambda1 ~ 3 pi / int_0^1 f1 dx (Eq. B8), f1(0.5).
4. Table 1 ratios and the e-folding time 1/lambda1.
5. The late-time amplitude c1 (sigma ~ c1 sigma0 f1 e^(-lambda1 t/t0), f1(1) = 1) by projecting the discrete
   uniform initial traction on the discrete modes (not in the key; the graders of MP-02, MP-05, MP-09 quote 1.531).
6. The round's automatic checker (E/checks.py, unchanged) on the saved results of the nine answered runs, and the
   count of '*_ok' flags that summary_data.py reports as 'key checks'.
7. MP-44's saved solve.py, replayed at N = 40 as delivered and with T_k(cos 2 theta) changed to T_k(-cos 2 theta).
Run from this folder:  python -I grain-boundary_check.py   (numpy, scipy)
"""
import datetime
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy
import scipy.linalg as sl
from scipy.integrate import quad

BUNDLE = Path(__file__).resolve().parents[2]
EV_MSB = BUNDLE / "<study-folder>/multiphysics_round"
EV_RUNS = BUNDLE / "evidence/corpus-runs/multiphysics"
pi = math.pi

print(f"# grain-boundary_check.py run {datetime.datetime.now(datetime.timezone.utc):%Y-%m-%dT%H:%M:%SZ} "
      f"with python {sys.version.split()[0]} numpy {np.__version__} scipy {scipy.__version__}")
print("# command: python -I grain-boundary_check.py (from assets/solution-evidence-2026-10-25)")

TABLE1 = [11.74265, 264.6348, 1357.440, 3845.215, 8328.149]       # Gao et al. 1999, Table 1, direct method


def S(x, xi):
    return 1 / (x - xi) - 1 / (x + xi) - 2 * xi * (x - xi) / (x + xi) ** 3


def Kker(x, xi):                                                   # Eq. 17: f(x) = int_0^1 K(x,xi) f'''(xi) dxi
    return np.where(x <= xi, 0.5 * xi ** 2 - 0.5 * (x - xi) ** 2, 0.5 * xi ** 2)


def grid(N):
    n = 2 * N + 1
    th = (2 * np.arange(1, N + 1) - 1) * pi / (2 * n)                # xi_j = cos(th_j), Eq. A4, positive half
    xi = np.cos(th)
    x = np.cos(np.arange(1, N + 1) * pi / n)                        # collocation points, Eq. A6
    return n, th, xi, x


def cheb_coeffs(th, r):
    """R(cos th) = sum_k a_k cos(2 k th), interpolating the N nodal values (R is even in xi)."""
    T = np.cos(2 * np.outer(th, np.arange(len(th))))
    return np.linalg.solve(T, r)


def R_of(a, theta):
    return np.cos(2 * np.outer(np.atleast_1d(theta), np.arange(len(a)))) @ a


def f_of(a, x):
    """f(x) = int_0^1 K(x,xi) R(xi)/sqrt(1-xi^2) dxi = int_0^{pi/2} K(x, cos th) R(cos th) dth."""
    g = lambda t: float(Kker(x, math.cos(t)) * R_of(a, t)[0])
    tx = math.acos(min(max(x, 0.0), 1.0))
    pts = [tx] if 0 < tx < pi / 2 else None
    return quad(g, 0, pi / 2, points=pts, limit=200, epsabs=1e-13, epsrel=1e-12)[0]


def minus_fpp(a, x):
    """-f''(x) = int_x^1 f''' dxi = int_0^{arccos x} R(cos th) dth: the shape of the wedge opening."""
    tx = math.acos(x)
    k = np.arange(1, len(a))
    return a[0] * tx + np.sum(a[1:] * np.sin(2 * k * tx) / (2 * k))


def eigen(N):
    n, th, xi, x = grid(N)
    X, XI = np.meshgrid(x, xi, indexing="ij")
    A = -(pi / n) * S(X, XI)                                        # Eq. A8
    B = (pi / n) * Kker(X, XI)
    w, V = sl.eig(A, B)
    ok = np.isfinite(w) & (np.abs(w.imag) < 1e-8 * np.abs(w.real)) & (w.real > 0)
    order = np.argsort(np.where(ok, w.real, np.inf))
    return n, th, xi, x, A, B, w, V, order


print("\n1. Eigenvalues by the paper's Gauss-Chebyshev collocation (Appendix A)")
store = {}
for N in (30, 60, 120):
    n, th, xi, x, A, B, w, V, order = eigen(N)
    lam = [w[order[i]].real for i in range(3)]
    r1 = V[:, order[0]].real
    a = cheb_coeffs(th, r1)
    f1, fh = f_of(a, 1.0), f_of(a, 0.5)
    store[N] = dict(lam=lam, a=a / f1, n=n, th=th, xi=xi, x=x, A=A, B=B, w=w, V=V, order=order)
    print(f"  N={N:3d}: lambda1 {lam[0]:.6f}  lambda2 {lam[1]:.4f}  lambda3 {lam[2]:.2f}  "
          f"f1(0.5)/f1(1) {fh / f1:.6f}")
lam1, lam2, lam3 = store[120]["lam"]
print(f"  Table 1 (paper): lambda1 {TABLE1[0]}  lambda2 {TABLE1[1]}  lambda3 {TABLE1[2]}")
print(f"  relative difference at N=120: lambda1 {lam1 / TABLE1[0] - 1:.1e}  lambda2 {lam2 / TABLE1[1] - 1:.1e}")

a120 = store[120]["a"]
print("  first mode at N=120, normalized f1(1) = 1:")
for xx in (0.1, 0.25, 0.5, 0.75, 0.9):
    print(f"    x={xx:4.2f}: f1 {f_of(a120, xx):.5f}   opening shape -f1'' {minus_fpp(a120, xx):.5f}   "
          f"crack shape 3 sqrt(1-x^2) {3 * math.sqrt(1 - xx * xx):.5f}")
print(f"    x=0.00: opening shape -f1'' {minus_fpp(a120, 0.0):.5f} (crack shape 3)")

print("\n2. Traction-free limit (edge crack of depth h), K/(sigma0 sqrt(pi h)) = pi |R(1)|")
for N in (30, 60, 120):
    n, th, xi, x = grid(N)
    X, XI = np.meshgrid(x, xi, indexing="ij")
    Rk = np.linalg.solve((pi / n) * S(X, XI), np.ones(N))
    a = cheb_coeffs(th, Rk)
    print(f"  N={N:3d}: K_inf/(sigma0 sqrt(pi h)) = {pi * abs(R_of(a, 0.0)[0]):.5f}")
print("  handbook edge-crack value quoted in the key: 1.1215 (not from the paper)")

print("\n3. Appendix B closed form (crack-shaped first mode)")
fB = lambda x: -0.25 * (6 * math.sqrt(1 - x * x) - 2 * (1 - x * x) ** 1.5 + 6 * x * math.asin(x) - 3 * pi * x - 4)
intB = quad(fB, 0, 1)[0]
print(f"  int_0^1 f1 dx = {intB:.6f}  (closed form 1 - 3 pi/32 = {1 - 3 * pi / 32:.6f})")
print(f"  Eq. B8: lambda1 ~ 3 pi / int f1 = {3 * pi / intB:.4f} (paper 13.3595); "
      f"{3 * pi / intB / TABLE1[0] - 1:+.1%} from Table 1")
print(f"  f1(0) = {fB(0):.3g}, f1(1) = {fB(1):.6f}, f1(0.5) = {fB(0.5):.4f}")
h = 1e-4
print(f"  check f1'' = -3 sqrt(1-x^2) at x=0.3: {(fB(0.3 + h) - 2 * fB(0.3) + fB(0.3 - h)) / h ** 2:.5f} "
      f"vs {-3 * math.sqrt(1 - 0.09):.5f}")

print("\n4. Decay ratios (Table 1) and the e-folding time")
print("  lambda_n/lambda_1 = " + ", ".join(f"{v / TABLE1[0]:.1f}" for v in TABLE1[1:]))
print(f"  1/lambda1 = {1 / TABLE1[0]:.5f} t0")

print("\n5. Late-time amplitude c1 of the first mode (f1(1) = 1), discrete projection of sigma(x,0) = 1")
for N in (30, 60, 120):
    d = store[N]
    n, xi, V, order, B = d["n"], d["xi"], d["V"], d["order"], d["B"]
    r0 = np.linalg.solve(B, np.ones(N))
    c = np.linalg.solve(V, r0.astype(complex))
    v1 = V[:, order[0]]
    f1_at_1 = np.sum((pi / n) * 0.5 * xi ** 2 * v1)
    print(f"  N={N:3d}: c1 = {(c[order[0]] * f1_at_1).real:.5f}")

print("\n6. The round's automatic checker on the saved results (E/checks.py, unchanged)")
summary = {r["label"]: r for r in json.loads((EV_MSB / "summary_data.json").read_text()) if r["problem"] == "E"}
for lab in ("MP-02", "MP-09", "MP-03", "MP-05", "MP-40", "MP-44", "MP-47", "MP-18", "MP-20"):
    out = subprocess.run([sys.executable, "-I", str(EV_MSB / "E/checks.py"), str(EV_RUNS / "runs" / lab)],
                         capture_output=True, text=True, check=True).stdout
    chk = json.loads(out)
    ok = [v for k, v in chk.items() if k.endswith("_ok") and v is not None]
    res = json.loads((EV_RUNS / "runs" / lab / "results.json").read_text())
    print(f"  {lab}: lambda1 {res['lambda1']:<20} class {chk['lambda1_class']:<11} '_ok' flags {sum(ok)}/{len(ok)} "
          f"(summary_data {summary[lab]['key_checks'][0]}/{summary[lab]['key_checks'][1]})  "
          f"own checks {chk['own_lambda2_gt_lambda1']}, {chk['own_time_consistent_with_lambda1']}")
print("  the '_ok' flags are lambda2, mode1_shape, K_inf and relaxation_time; lambda1_class is not among them")

print("\n7. MP-44's solve.py replayed at N = 40 (evidence copy, unchanged except where stated)")
th = 0.4
print(f"  identity: s = sin^2(th) gives 2s-1 = {2 * math.sin(th) ** 2 - 1:.6f} = -cos(2 th) = {-math.cos(2 * th):.6f}")
spec = importlib.util.spec_from_file_location("mp44", EV_RUNS / "runs/MP-44/solve.py")
mp44 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mp44)


def first_modes(N):
    lam, vec, Zc = mp44.solve(N)
    cand = sorted((l.real, i) for i, l in enumerate(lam) if abs(l.imag) < 1e-6 * max(1, abs(l.real)) and l.real > 1e-6)
    c = vec[:, cand[0][1]].real
    half = sum(c[j] * mp44.fk(j, 0.5) for j in range(N)) / sum(c[j] * mp44.fk(j, 1.0) for j in range(N))
    return cand[0][0], cand[1][0], half


l1, l2, hf = first_modes(40)
print(f"  as delivered, T_k(cos 2th):       lambda1 {l1:.5f}  lambda2 {l2:.3f}  f1(0.5) {hf:.4f}")
orig = mp44.theta_int


def theta_int_fixed(k, Z, hfun):
    thZ = np.arcsin(np.sqrt(Z)); t = 0.5 * thZ * (mp44.xth + 1.0); w = 0.5 * thZ * mp44.wth; s = np.sin(t) ** 2
    return 2.0 * np.sum(w * hfun(s) * mp44.chebT(k, -np.cos(2.0 * t)))


mp44.theta_int = theta_int_fixed
l1, l2, hf = first_modes(40)
print(f"  one sign changed, T_k(-cos 2th):  lambda1 {l1:.5f}  lambda2 {l2:.3f}  f1(0.5) {hf:.4f}")
mp44.theta_int = orig
