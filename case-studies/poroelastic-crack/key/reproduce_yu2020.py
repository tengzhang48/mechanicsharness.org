"""Reproduce Yu, Bouklas, Landis & Huang 2020 (J. Appl. Mech. 87, 031005), Eqs. 12-13, the
short-time limit of a stationary crack in an immersed linearly poroelastic gel (log P371).
Uses GPT-6.1-sol's own factorization (run MP-37, factor_check.py, unmodified; read-only copy
of the run folder) at nu_u = 0.5, eta = (nu_u - nu)/(1 - nu_u); Richardson step as in P370.
Yu et al.: K(0+) = K0 [1/(2(1-nu))]^alpha with alpha ~ 0.735; J*(0+) = J0 [1/(2(1-nu))]^(2 alpha - 1).
Run: python3 reproduce_yu2020.py > reproduce_yu2020.txt"""
import math, importlib.util
SRC = "<volume-2>/corpus-runs/multiphysics/runs/MP-37/factor_check.py"
spec = importlib.util.spec_from_file_location("fc", SRC); fc = importlib.util.module_from_spec(spec); spec.loader.exec_module(fc)
print("case            nu    nu_u  eta    K0+/K0    alpha   Yu Eq.12 (0.735)  J*/J0 here  Yu Eq.13")
for name, nu, nuu in [("B control", 0.2, 0.4)] + [("gel", n, 0.5) for n in (0.1, 0.2, 0.24, 0.3, 0.4, 0.45, 0.49)]:
    eta = (nuu - nu) / (1 - nuu)
    rows = [fc.factor(eta, n, 0.1) for n in (65536, 262144)]
    n11 = 2 * rows[-1]["n11"] - rows[-2]["n11"]
    K = n11 / (1 + eta); a = -math.log(K) / math.log(1 + eta)
    J = K ** 2 * (1 - nu) / (1 - nuu)
    yuK = (1 / (2 * (1 - nu))) ** 0.735 if nuu == 0.5 else float("nan")
    yuJ = (1 / (2 * (1 - nu))) ** (2 * 0.735 - 1) if nuu == 0.5 else float("nan")
    print(f"{name:12s} {nu:5.2f} {nuu:5.2f} {eta:6.3f} {K:9.6f} {a:7.4f} {yuK:12.6f} {J:14.4f} {yuJ:9.4f}", flush=True)
