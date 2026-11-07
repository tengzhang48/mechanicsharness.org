#!/usr/bin/env python3
"""Independent reviewer checks: ODE solve + energy quadrature, not submitted formula copy."""
import json
import math
import numpy as np
from scipy.integrate import solve_bvp, simpson

NU = 0.3
EBAR_OVER_E = 1.0 / (1.0 - NU**2)  # Ebar/E
# Work in dimensionless units: ell=1, Ebar=1, sigma0=1, h=1 => k=1

def strip_bvp(W):
    """Solve v'' - v = 0, v'(±W/2) = -1. Returns y, v, vp."""
    a, b = -W / 2.0, W / 2.0

    def fun(x, y):
        return np.vstack((y[1], y[0]))

    def bc(ya, yb):
        return np.array([ya[1] + 1.0, yb[1] + 1.0])

    n = max(80, int(400 * W))
    x = np.linspace(a, b, n)
    y0 = np.zeros((2, x.size))
    sol = solve_bvp(fun, bc, x, y0, tol=1e-10, max_nodes=8000)
    assert sol.success, sol.message
    return sol.x, sol.y[0], sol.y[1]


def U_from_fields(y, v, vp):
    """U/(h sigma0^2 ell / Ebar) with ell=h=Ebar=sigma0=1: int (vp + 0.5 vp^2 + 0.5 v^2) dy"""
    dens = vp + 0.5 * vp**2 + 0.5 * v**2
    return float(simpson(dens, x=y))


def gamma_from_U(U_ahead, U_behind, N=1.0):
    # G = (U_ahead - U_behind)/(N h); gamma = G Ebar/(sigma0^2 ell)
    return (U_ahead - U_behind) / N


rows_ref = []
ell_over_H = [0.5, 1.0, 1.5, 2.0]
independent = []
for r in ell_over_H:
    w = 1.0 / r  # H/ell
    yH, vH, vpH = strip_bvp(w)
    y2, v2, vp2 = strip_bvp(2 * w)
    UH = U_from_fields(yH, vH, vpH)
    U2 = U_from_fields(y2, v2, vp2)
    gA = gamma_from_U(0.0, UH, 1.0)
    gB = gamma_from_U(U2, 2 * UH, 1.0)
    kA = math.sqrt((1 - NU**2) * gA)
    kB = math.sqrt((1 - NU**2) * gB)
    independent.append(
        {
            "ell_over_H": float(r),
            "H_over_ell": float(w),
            "U_H": float(UH),
            "U_2H": float(U2),
            "gamma_A": float(gA),
            "gamma_B": float(gB),
            "kI_A": float(kA),
            "kI_B": float(kB),
        }
    )

# submitted
sub = json.load(open("submission/results.json"))
tol = 1e-6
agree = 0
details = []
for i, row in enumerate(sub["rows"]):
    ind = independent[i]
    for key in ("gamma_A", "gamma_B", "kI_A", "kI_B"):
        s = row[key] if key != "kI_A" else row["kI_A"]
        if key == "kI_A":
            s = row["kI_A"]
        elif key == "kI_B":
            s = row["kI_B"]
        else:
            s = row[key]
        d = abs(s - ind[key])
        ok = bool(float(d) < tol)
        if ok:
            agree += 1
        details.append({"ell_over_H": float(row["ell_over_H"]), "qty": key, "sub": float(s), "ind": float(ind[key]), "abs_diff": float(d), "agree": ok})

extra = {}
y, v, vp = strip_bvp(2.0)
d2 = np.gradient(np.gradient(v, y), y)
extra["max_ode_residual_fd"] = float(np.max(np.abs(d2 - v)))
extra["face_vp_plus"] = float(vp[-1] + 1)
extra["face_vp_minus"] = float(vp[0] + 1)
# independent tanh identity (not copied from submission file; derived from BVP energy)
# Also evaluate closed form from first principles: U = -tanh(W/2)
import math as _m
extra["closed_form_vs_bvp_max"] = max(
    abs(-_m.tanh(ind["H_over_ell"] / 2) - ind["U_H"]) for ind in independent
)

out = {
    "independent": independent,
    "agreement_count": agree,
    "agreement_of": 16,
    "tolerance": tol,
    "details": details,
    "extra": extra,
}
print(json.dumps(out, indent=2))
json.dump(out, open("review_independent_check_results.json", "w"), indent=2)
