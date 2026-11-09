"""Independent recomputation of the P038 T1 and T2 keys (audit 2026-10-25).

Does not import or copy the builder's reference.py.
T1: route A = paper Eqs. 45-46 (Bao & Suo 1992, p. 362, read from page image) with K_tip = K0
    (paper p. 359: "rising part of the R-curve is given by (17) by letting K_tip = K0";
    p. 361: "subjected to G_tip = Gamma0").
    route B = Abel/Paris double integral for the full gap of a crack under symmetric pressure,
    delta(x) = 4/(pi E') * int_{|x|}^{c} t/sqrt(t^2-x^2) * [int_{-t}^{t} p(s)/sqrt(t^2-s^2) ds] dt,
    which is not the log-kernel used by the builder. Inner integral done in closed form.
T2: own shooting (solve_ivp) of delta'''' = -C t(delta) on [0, L] with delta(0)=delta'(0)=0,
    delta''(0)=sqrt(2 C Gamma0), unknowns (L, delta'''(0)); plus continuation in the rear opening
    to check that L and M increase monotonically up to the event, and that delta is monotone on [0, L].
Also computes the plausible wrong answers used as mutants.
Run: <volume>/venvs/coupfe-eda/bin/python -I recompute_P038.py
"""
import json, math, os
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, fsolve

out = {}

# ---------------- T1 ----------------
E, nu, a, sb, dc, G0 = 130e9, 0.28, 2.40e-3, 160e6, 9.00e-6, 22.0
Ep = E / (1 - nu**2)
K0 = math.sqrt(Ep * G0)


def Sigma_tip(c, K0=K0, sb=sb, a=a):
    # Eq. 46 with K_tip = K0, solved for Sigma
    return K0 / math.sqrt(math.pi * c) + 2 * sb / math.pi * math.acos(a / c)


def delta_eq45(c, S, Ep=Ep, sb=sb, a=a):
    eta = c / a
    return sb * a / Ep * (4 * math.sqrt(eta**2 - 1) * (S / sb - 2 / math.pi * math.acos(1 / eta))
                          + 8 / math.pi * math.log(eta))


def solve_T1(Ep=Ep, sb=sb, a=a, dc=dc, K0=K0, opening_factor=1.0):
    f = lambda c: opening_factor * delta_eq45(c, Sigma_tip(c, K0, sb, a), Ep, sb, a) - dc
    c = brentq(f, a * (1 + 1e-12), a * 1e4, xtol=1e-18, rtol=1e-15)
    return Sigma_tip(c, K0, sb, a), c - a


S1, L1 = solve_T1()


# route B: Abel/Paris representation (independent of Eq. 45)
def inner(tt, S, c):
    # int_{-t}^{t} p(s)/sqrt(t^2-s^2) ds with p = S on |s|<a, S - sb on a<|s|<c (t <= c)
    if tt <= a:
        return S * math.pi
    return S * math.pi - sb * 2 * (math.pi / 2 - math.asin(a / tt))


def delta_abel(x, c, S):
    g = lambda tt: tt / math.sqrt(tt**2 - x**2) * inner(tt, S, c)
    # substitution tt = x cosh(w) removes the endpoint singularity
    wmax = math.acosh(c / x)
    h = lambda w: inner(x * math.cosh(w), S, c) * x * math.cosh(w)
    val, err = quad(h, 0.0, wmax, epsabs=0, epsrel=1e-13, limit=200)
    return 4 / (math.pi * Ep) * val


def K_weight(c, S):
    val, _ = quad(lambda s: (S if s < a else S - sb) / math.sqrt(c**2 - s**2), 0, c,
                  points=[a], epsabs=0, epsrel=1e-13, limit=200)
    return 2 * math.sqrt(c / math.pi) * val


def solve_T1_B():
    def F(c):
        # Sigma from integrated tip condition (linear in S): K(S) = K(0) + S*sqrt(pi c)
        Kb = K_weight(c, 0.0)
        S = (K0 - Kb) / math.sqrt(math.pi * c)
        return delta_abel(a, c, S) - dc, S
    c = brentq(lambda c: F(c)[0], a * 1.01, a * 10, xtol=1e-16, rtol=1e-14)
    return F(c)[1], c - a


S1b, L1b = solve_T1_B()
# branch monotonicity: Sigma increasing on (a, c*]
cs = np.linspace(a * (1 + 1e-6), a + L1, 4001)
Ss = np.array([Sigma_tip(c) for c in cs])
ds_open = np.array([delta_eq45(c, Sigma_tip(c)) for c in cs])
# load maximum on the full branch (if any) for context
cfine = a * np.exp(np.linspace(1e-6, 6, 20001))
Sf = np.array([Sigma_tip(c) for c in cfine])
out["T1"] = {
    "Sigma_MPa_eq45_46": S1 / 1e6, "L_mm_eq45_46": L1 * 1e3,
    "Sigma_MPa_abel_route": S1b / 1e6, "L_mm_abel_route": L1b * 1e3,
    "rel_diff_Sigma": abs(S1b - S1) / S1, "rel_diff_L": abs(L1b - L1) / L1,
    "key_Sigma_MPa": 133.730286053811, "key_L_mm": 4.08560913824553,
    "rel_err_vs_key_Sigma": abs(S1 / 1e6 - 133.730286053811) / 133.730286053811,
    "rel_err_vs_key_L": abs(L1 * 1e3 - 4.08560913824553) / 4.08560913824553,
    "Sigma_monotone_to_event": bool(np.all(np.diff(Ss) > 0)),
    "root_opening_monotone_to_event": bool(np.all(np.diff(ds_open) > 0)),
    "initiation_Sigma_MPa": Sigma_tip(a * (1 + 1e-15)) / 1e6,
    "branch_max_Sigma_MPa_up_to_c_over_a_e6": float(Sf.max() / 1e6),
    "branch_max_at_L_mm": float((cfine[Sf.argmax()] - a) * 1e3),
    "a_over_material_length": a / (Ep * dc / sb),
}

# T1 mutants (plausible wrong answers)
mut = {}
S, L = solve_T1(K0=0.0); mut["Gamma0_omitted"] = (S / 1e6, L * 1e3)
S, L = solve_T1(opening_factor=0.5); mut["half_gap_as_full_gap(one-face disp = delta_c)"] = (S / 1e6, L * 1e3)
mut["total_L_reported(2L)"] = (S1 / 1e6, 2 * L1 * 1e3)
S, L = solve_T1(Ep=E, K0=math.sqrt(E * G0)); mut["plane_stress_E"] = (S / 1e6, L * 1e3)
S, L = solve_T1(K0=math.sqrt(Ep * 2 * G0)); mut["Gamma0_doubled"] = (S / 1e6, L * 1e3)
mut["SSB_LEFM_plateau_no_zone"] = (math.sqrt(Ep * (G0 + sb * dc) / (math.pi * a)) / 1e6, None)
mut["seen_R03_without_dsaz"] = (158.524900100573, 34.324303918593)
out["T1_mutants"] = mut

# ---------------- T2 ----------------
b, E1, h1, E2, h2 = 25.0e-3, 110e9, 0.60e-3, 190e9, 0.90e-3
sb2, dc2, G02 = 0.20e6, 0.20e-3, 1.40


def T2_solve(D1, D2, G0=G02, sb=sb2, dc=dc2, rear=1.0, guess=(1.68, 1.18)):
    C = 1 / D1 + 1 / D2
    ell = (dc / (C * sb)) ** 0.25
    kap = math.sqrt(2 * C * G0) * ell**2 / dc  # dimensionless u''(0)

    def shoot(lam, A):
        sol = solve_ivp(lambda z, y: [y[1], y[2], y[3], y[0] - 1.0], (0, lam), [0, 0, kap, A],
                        method="DOP853", rtol=1e-13, atol=1e-15, dense_output=True)
        return sol

    def res(p):
        lam, A = p
        y = shoot(lam, A).y[:, -1]
        return [y[0] - rear, y[3]]

    lam, A = fsolve(res, guess, xtol=1e-14)
    sol = shoot(lam, A)
    yL = sol.y[:, -1]
    m = yL[2] * dc / ell**2 / C          # moment per width from rear curvature
    zz = np.linspace(0, lam, 2001)
    U = sol.sol(zz)
    return {"lam": lam, "A": A, "L_mm": lam * ell * 1e3, "M_Nm": m * b,
            "M_energy_Nm": b * math.sqrt(2 * (G0 + sb * dc * (rear - rear**2 / 2)) / C),
            "dmid_um": sol.sol(lam / 2)[0] * dc * 1e6, "C": C, "ell_mm": ell * 1e3,
            "res": res((lam, A)), "u_monotone": bool(np.all(np.diff(U[0]) >= -1e-14)),
            "u_min": float(U[0].min()), "u_max": float(U[0].max())}


D1, D2 = E1 * h1**3 / 12, E2 * h2**3 / 12
r = T2_solve(D1, D2)
# continuation in rear opening: L, M must increase up to rear = 1
cont = []
g = (0.3, 0.4)
for ur in np.linspace(0.02, 1.0, 50):
    rr = T2_solve(D1, D2, rear=ur, guess=g)
    g = (rr["lam"], rr["A"])
    cont.append((ur, rr["L_mm"], rr["M_Nm"], rr["u_monotone"]))
Ls = [c[1] for c in cont]; Ms = [c[2] for c in cont]
out["T2"] = {"M_Nm_shoot": r["M_Nm"], "M_Nm_energy": r["M_energy_Nm"], "L_mm": r["L_mm"],
             "delta_mid_um": r["dmid_um"], "lambda": r["lam"], "residual": list(map(float, r["res"])),
             "profile_monotone": r["u_monotone"],
             "key": {"M": 0.212625787376721, "L": 10.8001023708973, "dmid": 38.2188979260933},
             "rel_err_vs_key": {"M": abs(r["M_Nm"] - 0.212625787376721) / 0.212625787376721,
                                "L": abs(r["L_mm"] - 10.8001023708973) / 10.8001023708973,
                                "dmid": abs(r["dmid_um"] - 38.2188979260933) / 38.2188979260933},
             "continuation_L_increasing": bool(np.all(np.diff(Ls) > 0)),
             "continuation_M_increasing": bool(np.all(np.diff(Ms) > 0)),
             "continuation_profiles_monotone": bool(all(c[3] for c in cont)),
             "first_continuation_point": cont[0][:3]}

m2 = {}
rr = T2_solve(D1, D2, G0=0.0, guess=(1.87, 1.0)); m2["Gamma0_omitted"] = (rr["M_Nm"], rr["L_mm"], rr["dmid_um"])
Dav = ((E1 + E2) / 2) * ((h1 + h2) / 2) ** 3 / 12
rr = T2_solve(Dav, Dav); m2["averaged_equal_arms"] = (rr["M_Nm"], rr["L_mm"], rr["dmid_um"])
rr = T2_solve(D1, 1e30); m2["only_upper_arm_compliance"] = (rr["M_Nm"], rr["L_mm"], rr["dmid_um"])
m2["moment_per_width_reported_as_M(N)"] = (r["M_Nm"] / b, r["L_mm"], r["dmid_um"])
m2["sum_of_two_moments"] = (2 * r["M_Nm"], r["L_mm"], r["dmid_um"])
m2["delta_half_gap(one arm share)"] = (r["M_Nm"], r["L_mm"], r["dmid_um"] / 2)
out["T2_mutants"] = m2

# ---------------- T3 numbers (optional diagnostics) ----------------
EpA = 120e9 / (1 - 0.25**2)
out["T3"] = {"A_W": 0.5 * 240e6 * 1.5e-6, "A_lm_mm": EpA * 1.5e-6 / 240e6 * 1e3,
             "A_Kc": math.sqrt(EpA * 180) / 1e6, "B_W": 0.5 * 300e6 * 30e-6,
             "C_Wfull": 2 / 3 * 8e6 * 80e-6, "C_Wcap": 2 / 3 * 8e6 * 80e-6 * (3 / 80) ** 1.5,
             "C_frac": (3 / 80) ** 1.5, "C_over_G0": 2 / 3 * 8e6 * 80e-6 * (3 / 80) ** 1.5 / 12}

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "recompute_P038_results.json"), "w") as fh:
    json.dump(out, fh, indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
