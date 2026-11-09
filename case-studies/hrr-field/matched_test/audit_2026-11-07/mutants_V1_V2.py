"""Wrong answers for the checker test (audit 2026-11-07): the errors seen in runs, plus plausible ones.

V1 mutants are computed with this folder's independent field (recompute_V1.py): same eigenfield,
wrong amplitude/convention.  V2 mutants use closed forms of the leading (T, nu_t) field, each with one
deliberate error.  Writes mutants.json (consumed by check_V1_V2.py --audit) and mutants.out.
"""
import json
import os
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import recompute_V1 as v1                       # noqa: E402  (module-level shooting only)
from recompute_V2 import full_law               # noqa: E402

M = {}
lines = []


def add(name, prob, vals, note):
    M[name] = {"problem": prob, "values": {k: float(v) for k, v in vals.items()}, "note": note}
    lines.append(f"{prob} {name}: {note}\n    " + ", ".join(f"{k}={float(v):.6g}" for k, v in vals.items()))


# ---------------------------------------------------------------- V1 ---------------------------
n, c, mu, lam, J = v1.n, v1.c, v1.mu, v1.lam, v1.J
J1 = v1.J_circle(1.0, r=1.0)
A = (J / J1) ** (1 / (n + 1))
base = v1.outputs(A)


def polar_J(A_, wrong_parity=False, r=1.0):
    """Circle J with the field split into polar parts; wrong_parity flips gamma_theta in the lower
    half inside w_,x only (the M20 error: f' extended odd instead of even)."""
    def f(th):
        g, gp, h = v1.ang(th)
        tr, tt = A_ * gp, -mu * A_ * g                     # r^(mu-1) = 1 at r = 1
        tm = np.hypot(tr, tt)
        gr_, gt_ = c * tm ** (n - 1) * tr, c * tm ** (n - 1) * tt
        if wrong_parity and th < 0:
            gt_ = -gt_
        W = n / (n + 1) * tm * c * tm ** n
        wx = gr_ * np.cos(th) - gt_ * np.sin(th)
        return r * (W * np.cos(th) - tr * wx)
    return quad(f, -np.pi, 0, limit=200, epsrel=1e-12)[0] + quad(f, 0, np.pi, limit=200, epsrel=1e-12)[0]


Jp = polar_J(1.0, wrong_parity=True)
lines.append(f"V1 parity check: J(A=1) correct {polar_J(1.0):.10e}, wrong parity {Jp:.10e}, "
             f"ratio {Jp / polar_J(1.0):.6f} (M20 grade: 624.1363418/953.8687488 = {624.1363418/953.8687488:.6f})")
add("V1_parity_error_lower_half_fprime_odd", "V1", v1.outputs((J / Jp) ** (1 / (n + 1))),
    "f' extended with the wrong parity in J (computed); M20's error class")
for drop in ("traction", "W"):
    Jd = v1.J_circle(1.0, r=1.0, drop=drop)
    if Jd > 0:
        add(f"V1_J_term_dropped_{drop}", "V1", v1.outputs((J / Jd) ** (1 / (n + 1))),
            f"J with the {'t_z w_,x' if drop == 'traction' else 'W N_x'} term dropped (J1 ratio {Jd / J1:.4f})")
    else:
        lines.append(f"V1 J_term_dropped_{drop}: J(A=1) = {Jd:.3e} <= 0, no real amplitude (not formable)")
Jc = v1.J_circle(1.0, r=1.0, Wfun=lambda g: v1.tau_of_gamma(g) * g / (n + 1))
add("V1_complementary_energy_in_J", "V1", v1.outputs((J / Jc) ** (1 / (n + 1))),
    f"W replaced by complementary work tau*gamma/(n+1) (J1 ratio {Jc / J1:.4f})")
Jl = v1.J_circle(1.0, r=1.0, Wfun=lambda g: 0.5 * v1.tau_of_gamma(g) * g)
add("V1_linear_W_half_tau_gamma", "V1", v1.outputs((J / Jl) ** (1 / (n + 1))),
    f"W = tau.gamma/2 (linear-elastic habit) (J1 ratio {Jl / J1:.4f})")
add("V1_law_applied_to_tensor_shear", "V1", v1.outputs((J / (2 * J1)) ** (1 / (n + 1)), cfac=2.0),
    "law read as eps = (gamma0/tau0)(tau/tau0)^(n-1) tau with eps = gamma/2: compliance doubled")
rc = {}
for th, k in ((0.0, "r_0"), (np.pi / 2, "r_90"), (np.pi, "r_180")):
    g, gp, h = v1.ang(th); Q = np.hypot(gp, mu * g)
    rc[k] = ((0.02 / c) ** (1 / n) / (A * Q)) ** (1 / (mu - 1))
add("V1_contour_on_tensor_shear", "V1", rc, "contour |eps_tensor| = 0.01, i.e. |gamma| = 0.02")
add("V1_tau_thetaz_sign_flipped", "V1", {"tau_tz": -base["tau_tz"]}, "sign of tau_theta z reversed")
add("V1_tau_rz_sign_flipped", "V1", {"tau_rz": -base["tau_rz"]}, "sign of tau_rz reversed")
t60 = np.pi / 3
add("V1_cartesian_components_reported", "V1",
    {"tau_rz": base["tau_rz"] * np.cos(t60) - base["tau_tz"] * np.sin(t60),
     "tau_tz": base["tau_rz"] * np.sin(t60) + base["tau_tz"] * np.cos(t60)},
    "tau_xz, tau_yz reported instead of polar components")
add("V1_slip_one_face_only", "V1", {"slip_0.02": base["slip_0.02"] / 2, "slip_0.08": base["slip_0.08"] / 2},
    "w(r, pi) reported instead of w(r, pi) - w(r, -pi)")
G = v1.tau0 / v1.gam0
K3 = np.sqrt(2 * G * J)
lin = {"slip_0.02": 4 * K3 / G * np.sqrt(0.02 / (2 * np.pi)), "slip_0.08": 4 * K3 / G * np.sqrt(0.08 / (2 * np.pi)),
       "tau_rz": K3 / np.sqrt(2 * np.pi * 0.02) * np.sin(t60 / 2), "tau_tz": K3 / np.sqrt(2 * np.pi * 0.02) * np.cos(t60 / 2)}
rl = (K3 / (G * 0.01)) ** 2 / (2 * np.pi)
lin.update({"r_0": rl, "r_90": rl, "r_180": rl})
add("V1_wrong_branch_linear_KIII_G_tau0_over_gamma0", "V1", lin,
    "linear mode III field (n = 1 branch) with G = tau0/gamma0, J = K^2/(2G)")
add("V1_M20_as_reported", "V1", {"slip_0.02": 1.83e-3, "slip_0.08": 2.49e-3, "tau_rz": 42.4, "tau_tz": 151.7,
                                 "r_0": 0.07167, "r_90": 0.03831, "r_180": 0.02048},
    "M20 ANSWER.md:15,16,20,30-32")

# ---------------------------------------------------------------- V2 ---------------------------
E, nu, Jv, Czz, sy = 70000.0, 0.30, 0.800, 0.005, 10.0
nut = lambda T: 0.5 - (0.5 - nu) * T / E
cs, sn = np.cos(np.pi / 6), np.sin(np.pi / 6)
ang60 = np.array([cs * (1 + sn * sn), cs ** 3, sn * cs * cs])


def v2_all(T, nu_, K, plane_strain=False, kappa=None):
    S = K / np.sqrt(2 * np.pi) * ang60
    if plane_strain:
        Crr = ((1 - nu_ ** 2) * S[0] - nu_ * (1 + nu_) * S[1]) / T
        Ctt = ((1 - nu_ ** 2) * S[1] - nu_ * (1 + nu_) * S[0]) / T
        kap = 3 - 4 * nu_
        pr = K / (np.sqrt(2 * np.pi) * T) * np.array([(1 + nu_) * (1 - 2 * nu_)] * 2 + [0.0])
        se = abs(1 - 2 * nu_)
    else:
        Crr, Ctt = (S[0] - nu_ * S[1]) / T, (S[1] - nu_ * S[0]) / T
        kap = (3 - nu_) / (1 + nu_)
        pr = K / (np.sqrt(2 * np.pi) * T) * np.array([1 - nu_, 1 - nu_, -2 * nu_])
        se = 1.0
    if kappa is not None:
        kap = kappa
    Cg = 2 * (1 + nu_) * S[2] / T
    Cy = K * (1 + nu_) * kap / (2 * np.sqrt(np.pi) * T)
    rmin = (np.max(np.abs(pr)) / 0.05) ** 2
    rmax = min((se * K) ** 2 / (2 * np.pi * (10 * sy) ** 2), 0.200)
    return {"T": T, "C_rr": Crr, "C_tt": Ctt, "C_g": Cg, "C_y": Cy, "r_min": rmin, "r_max": rmax}


T0 = brentq(lambda T: 2 * nut(T) * np.sqrt(Jv / (2 * np.pi * T)) - Czz, 1.0, E, xtol=1e-12)
K0 = np.sqrt(Jv * T0)
ref = v2_all(T0, nut(T0), K0)
lines.append("V2 closed-form base (should equal key): " + ", ".join(f"{k}={v:.10g}" for k, v in ref.items()))

Ta = brentq(lambda T: 2 * nu * np.sqrt(Jv / (2 * np.pi * T)) - Czz, 1.0, E)
add("V2_nu_held_at_0.30_E_replaced_by_T", "V2", v2_all(Ta, nu, np.sqrt(Jv * Ta)), "builder's named pitfall")
Tb = Jv / (2 * np.pi * Czz ** 2)
add("V2_incompressible_nu_t_half_paper_eq48", "V2", v2_all(Tb, 0.5, np.sqrt(Jv * Tb)),
    "nu_t = 1/2, i.e. eps_zz = -sigma/T as in the paper's Eq. 48 (elastic volume term dropped)")
Tc = brentq(lambda T: 2 * nut(T) * np.sqrt(Jv * T / (1 - nut(T) ** 2)) / (np.sqrt(2 * np.pi) * T) - Czz, 1.0, E)
add("V2_plane_strain_J_K_relation", "V2", v2_all(Tc, nut(Tc), np.sqrt(Jv * Tc / (1 - nut(Tc) ** 2))),
    "J = (1 - nu_t^2) K^2 / T (plane strain) with the plane-stress C_zz")
pe = v2_all(T0, nut(T0), K0, plane_strain=True)
pe.pop("T")
add("V2_plane_strain_field_T_correct", "V2", pe,
    "T right, then plane-strain compliance, kappa = 3 - 4 nu_t, sigma_zz = 2 nu_t sigma in the screens "
    "(interval empty: r_max < r_min)")
Te = brentq(lambda T: 2 * nut(T) * np.sqrt(Jv * E) / (np.sqrt(2 * np.pi) * T) - Czz, 1.0, E)
add("V2_elastic_E_in_J_K_relation", "V2", v2_all(Te, nut(Te), np.sqrt(Jv * E)), "J = K^2/E (elastic E, not T)")
add("V2_tensor_shear_reported", "V2", {"C_g": ref["C_g"] / 2}, "eps_r_theta coefficient reported instead of 2 eps_r_theta")
add("V2_C_g_sign_flipped", "V2", {"C_g": -ref["C_g"]}, "shear sign reversed")
add("V2_C_rr_C_tt_swapped", "V2", {"C_rr": ref["C_tt"], "C_tt": ref["C_rr"]}, "polar components swapped")
add("V2_C_y_plane_strain_kappa_only", "V2", {"C_y": v2_all(T0, nut(T0), K0, kappa=3 - 4 * nut(T0))["C_y"]},
    "kappa = 3 - 4 nu_t in the displacement only")
add("V2_lower_screen_trace_term_dropped_M08", "V2", {"r_min": ((1 + nut(T0)) * K0 / (np.sqrt(2 * np.pi) * T0) / 0.05) ** 2},
    "in-plane strains (1 + nu_t) sigma / T (M08's error)")
add("V2_lower_screen_in_plane_only", "V2", {"r_min": ((1 - nut(T0)) * K0 / (np.sqrt(2 * np.pi) * T0) / 0.05) ** 2},
    "eps_zz left out of the principal strains")
fl = lambda r: np.max(np.abs(np.linalg.eigvalsh(full_law(K0 / np.sqrt(2 * np.pi * r) * np.diag([1.0, 1.0, 0.0]), T0))))
add("V2_lower_screen_full_law_at_leading_stress", "V2", {"r_min": brentq(lambda r: fl(r) - 0.05, 1e-5, 1.0, xtol=1e-15)},
    "complete law (with the sigma_y term) evaluated at leading-field stresses (alternative reading; M22 ANSWER.md:79-80 mentions it)")
add("V2_upper_screen_sigma_y_not_10sigma_y", "V2", {"r_max": min(K0 ** 2 / (2 * np.pi * sy ** 2), 0.2)},
    "sigma_e >= sigma_y (factor 10 dropped): the 0.200 mm cap binds")
add("V2_M08_as_reported", "V2", {"T": 4816.47, "C_rr": 0.003942, "C_tt": 0.000633, "C_g": 0.005731,
                                 "C_y": 0.009139, "r_min": 0.02336, "r_max": 0.06133}, "M08 ANSWER.md:9-14")

json.dump(M, open(os.path.join(HERE, "mutants.json"), "w"), indent=1)
open(os.path.join(HERE, "mutants.out"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
