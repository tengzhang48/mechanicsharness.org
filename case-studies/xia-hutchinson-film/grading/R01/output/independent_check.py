#!/usr/bin/env python3
"""Independent 1D wake energy check (not solver FEM, not reference reproduce.py)."""
import json, math, os

nu = 0.3
one_m_nu2 = 1 - nu * nu
samples = [0.5, 1.0, 1.5, 2.0]


def gamma_A(ell_over_H):
    H_over_ell = 1.0 / ell_over_H
    return math.tanh(H_over_ell / 2.0)


def gamma_B(ell_over_H):
    H_over_ell = 1.0 / ell_over_H
    return 2.0 * math.tanh(H_over_ell / 2.0) - math.tanh(H_over_ell)


def kI(gamma):
    return math.sqrt(one_m_nu2 * gamma)


# Finite-difference energy of 1D ODE: v''=v, v'(0)=v'(L)=-1
def fd_gamma(L, n=2000):
    h = L / n
    # second-order FD with Neumann at both ends via ghost points
    # (v_{i+1}-2v_i+v_{i-1})/h^2 = v_i
    # left: (v1 - v_{-1})/(2h) = -1 => v_{-1} = v1 + 2h
    # right: (v_{n+1}-v_{n-1})/(2h) = -1 => v_{n+1} = v_{n-1} - 2h
    import numpy as np
    A = np.zeros((n + 1, n + 1))
    b = np.zeros(n + 1)
    for i in range(1, n):
        A[i, i - 1] = 1.0 / h**2
        A[i, i] = -2.0 / h**2 - 1.0
        A[i, i + 1] = 1.0 / h**2
    # i=0
    A[0, 0] = -2.0 / h**2 - 1.0
    A[0, 1] = 2.0 / h**2
    b[0] = -2.0 / h  # from ghost v_{-1}=v1+2h
    # i=n
    A[n, n] = -2.0 / h**2 - 1.0
    A[n, n - 1] = 2.0 / h**2
    b[n] = 2.0 / h  # ghost v_{n+1}=v_{n-1}-2h
    v = np.linalg.solve(A, b)
    # energy I = int [v' + (v'^2+v^2)/2] d eta
    # Simpson-like trapezoid of density; use centered v'
    vp = np.zeros_like(v)
    vp[0] = (v[1] - (v[1] + 2 * h)) / (2 * h)  # should be -1
    vp[-1] = ((v[-2] - 2 * h) - v[-2]) / (2 * h)
    vp[1:-1] = (v[2:] - v[:-2]) / (2 * h)
    dens = vp + 0.5 * (vp**2 + v**2)
    I = np.trapezoid(dens, dx=h)
    # G h / (h sigma0^2 ell / Ebar) = -I  for spacing S=L ell, one cell
    return -I


rows = []
for r in samples:
    gA = gamma_A(r)
    gB = gamma_B(r)
    L_H = 1.0 / r  # H/ell
    gA_fd = fd_gamma(L_H)
    gB_fd = 2 * fd_gamma(L_H) - fd_gamma(2 * L_H)
    rows.append(
        {
            "ell_over_H": r,
            "gamma_A_exact": gA,
            "gamma_B_exact": gB,
            "kI_A_exact": kI(gA),
            "kI_B_exact": kI(gB),
            "gamma_A_fd": gA_fd,
            "gamma_B_fd": gB_fd,
            "fd_err_A": abs(gA_fd - gA),
            "fd_err_B": abs(gB_fd - gB),
        }
    )

# small-H series
series = {
    "gamma_A_small": 0.5 / 2.0,  # H/ell=0.5 => gamma~ H/(2ell)=0.25 wait
}
# H/ell -> 0: gamma_A ~ H/(2 ell), at ell/H=2 => H/ell=0.5, ~0.25
# gamma_B ~ (H/ell)^3/4 = 0.5**3/4 = 0.03125

submitted_path = os.path.join(os.path.dirname(__file__), "submission", "results.json")
with open(submitted_path) as f:
    sub = json.load(f)

tol = 1e-3  # scientific table agreement (not 1e-6; FEM vs analytic)
n_agree = 0
cmp = []
for srow, erow in zip(sub["rows"], rows):
    for key, ek in [
        ("gamma_A", "gamma_A_exact"),
        ("gamma_B", "gamma_B_exact"),
        ("kI_A", "kI_A_exact"),
        ("kI_B", "kI_B_exact"),
    ]:
        sv = srow[key]
        ev = erow[ek]
        d = abs(sv - ev)
        ok = d <= 1e-3
        if ok:
            n_agree += 1
        cmp.append({"ell_over_H": srow["ell_over_H"], "q": key, "sub": sv, "exact": ev, "absdiff": d, "agree_1e-3": ok})

out = {
    "ell": "sqrt(Ebar h / k)",
    "gamma_A": "tanh(H/(2 ell))",
    "gamma_B": "2 tanh(H/(2 ell)) - tanh(H/ell)",
    "kI": "sqrt((1-nu^2) gamma)",
    "isolated_kI": math.sqrt(one_m_nu2),
    "rows": rows,
    "small_H": {
        "ell_over_H_2_gamma_A_series": 0.5 / 2.0,
        "ell_over_H_2_gamma_B_series": (0.5**3) / 4.0,
        "exact_A": gamma_A(2.0),
        "exact_B": gamma_B(2.0),
    },
    "comparison_submitted": cmp,
    "n_agree_1e-3": n_agree,
    "n_agree_1e-6": sum(1 for c in cmp if c["absdiff"] <= 1e-6),
    "n_compared": 16,
}
print(json.dumps(out, indent=2))
with open(os.path.join(os.path.dirname(__file__), "independent_check_output.json"), "w") as f:
    json.dump(out, f, indent=2)
print("n_agree_1e-3", n_agree, "n_agree_1e-6", out["n_agree_1e-6"])
