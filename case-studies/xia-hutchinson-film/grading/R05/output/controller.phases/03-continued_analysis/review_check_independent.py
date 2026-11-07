#!/usr/bin/env python3
"""Reviewer-owned deciding-quantity check.

Does not import submission/ or reference/ code. Discretizes the PROBLEM.md
incremental energy on a full-width 1D strip with natural total-traction BCs.

Dimensionless: v = Ebar u2 / (sigma0 ell), eta = x2/ell, L = S/ell.
I[v] = v(L)-v(0) + int_0^L (1/2)(v'^2 + v^2) d eta
Stationarity => v''=v and v'(0)=v'(L)=-1.
gamma(S) = -I for a spacing-S wake. Case A: gamma_A = -I(H).
Case B: gamma_B = I(2H) - 2 I(H).
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NU = 0.3
ONE_MINUS_NU2 = 1.0 - NU * NU
SAMPLES = [0.5, 1.0, 1.5, 2.0]
OUT_TXT = os.path.join(HERE, "review_check_independent_output.txt")
OUT_JSON = os.path.join(HERE, "review_check_independent_results.json")
SUB_JSON = os.path.join(HERE, "submission", "results.json")
REF_JSON = os.path.join(HERE, "reference", "reference.json")


def assemble_p1(L, n_el):
    n = n_el + 1
    h = L / n_el
    K = np.zeros((n, n))
    ke = (1.0 / h) * np.array([[1.0, -1.0], [-1.0, 1.0]]) + (
        h / 6.0
    ) * np.array([[2.0, 1.0], [1.0, 2.0]])
    for e in range(n_el):
        K[e : e + 2, e : e + 2] += ke
    rhs = np.zeros(n)
    rhs[0] = 1.0
    rhs[-1] = -1.0
    return K, rhs, h


def solve_p1(L, n_el):
    K, rhs, h = assemble_p1(L, n_el)
    v = np.linalg.solve(K, rhs)
    I = (v[-1] - v[0]) + 0.5 * float(v @ (K @ v))
    vp0 = (v[1] - v[0]) / h
    vpL = (v[-1] - v[-2]) / h
    return v, I, h, vp0, vpL


def linear_term_from_v(v, L):
    """int v' d eta = v(L)-v(0). Not identically zero on a symmetric strip."""
    return float(v[-1] - v[0])


def control_omit_foundation(L, n_el):
    Kfull, rhs, _ = assemble_p1(L, n_el)
    v = np.linalg.solve(Kfull, rhs)
    n = n_el + 1
    h = L / n_el
    Kdiff = np.zeros((n, n))
    ke = (1.0 / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])
    for e in range(n_el):
        Kdiff[e : e + 2, e : e + 2] += ke
    I = (v[-1] - v[0]) + 0.5 * float(v @ (Kdiff @ v))
    return -I


def control_incremental_traction(L, n_el):
    K, _, _ = assemble_p1(L, n_el)
    v = np.linalg.solve(K, np.zeros(n_el + 1))
    I = (v[-1] - v[0]) + 0.5 * float(v @ (K @ v))
    return -I, float(np.max(np.abs(v)))


def main():
    lines = []

    def log(s=""):
        print(s)
        lines.append(s)

    with open(SUB_JSON) as f:
        sub = json.load(f)
    with open(REF_JSON) as f:
        ref = json.load(f)
    sub_rows = {r["ell_over_H"]: r for r in sub["rows"]}
    ref_rows = {r["ell_over_H"]: r for r in ref["rows"]}

    log("Independent reviewer P1 FEM of PROBLEM.md strip energy")
    log("Natural BCs v'(0)=v'(L)=-1 (total traction). No submission functions.")
    log("nu=%.1f  1-nu^2=%.2f" % (NU, ONE_MINUS_NU2))

    log("")
    log("[CONV] L=1 (ell/H=1 Case A) vs tanh(0.5)")
    target = math.tanh(0.5)
    conv = []
    for n_el in (40, 80, 160, 320, 640, 1280):
        _, I, h, vp0, vpL = solve_p1(1.0, n_el)
        g = -I
        err = abs(g - target)
        conv.append(err)
        log(
            "  n_el=%4d  gamma=%.12f  |err|=%.3e  vp0=%.6f vpL=%.6f"
            % (n_el, g, err, vp0, vpL)
        )
    ratio = conv[-2] / conv[-1] if conv[-1] else float("inf")
    log("  err(640)/err(1280)=%.3f (expect ~4 for O(h^2))" % ratio)

    n_el = 4000
    log("")
    log("[MAIN] n_el=%d full-width P1 vs submission vs reference" % n_el)
    comparisons = []
    n_agree_1e12 = 0
    n_agree_1e6 = 0
    for eoh in SAMPLES:
        HoL = 1.0 / eoh
        vH, I_H, _, vp0, vpL = solve_p1(HoL, n_el)
        _, I_2H, _, _, _ = solve_p1(2.0 * HoL, n_el)
        gA = -I_H
        gB = -2.0 * I_H + I_2H
        lin = linear_term_from_v(vH, HoL)
        kA = math.sqrt(ONE_MINUS_NU2 * gA)
        kB = math.sqrt(ONE_MINUS_NU2 * gB)
        sA = sub_rows[eoh]
        rA = ref_rows[eoh]
        dA_sub = abs(gA - sA["gamma_A"])
        dB_sub = abs(gB - sA["gamma_B"])
        dKA_sub = abs(kA - sA["kI_A"])
        dKB_sub = abs(kB - sA["kI_B"])
        dA_ref = abs(gA - rA["gamma_A"])
        dB_ref = abs(gB - rA["gamma_B"])
        for d in (dA_sub, dB_sub, dKA_sub, dKB_sub):
            if d < 1e-12:
                n_agree_1e12 += 1
            if d < 1e-6:
                n_agree_1e6 += 1
        row = {
            "ell_over_H": eoh,
            "H_over_ell": HoL,
            "gamma_A_fem": gA,
            "gamma_B_fem": gB,
            "kI_A_fem": kA,
            "kI_B_fem": kB,
            "linear_term_I_vprime": lin,
            "vp0": vp0,
            "vpL": vpL,
            "d_gamma_A_sub": dA_sub,
            "d_gamma_B_sub": dB_sub,
            "d_kI_A_sub": dKA_sub,
            "d_kI_B_sub": dKB_sub,
            "d_gamma_A_ref": dA_ref,
            "d_gamma_B_ref": dB_ref,
            "sub_gamma_A": sA["gamma_A"],
            "sub_gamma_B": sA["gamma_B"],
            "sub_kI_A": sA["kI_A"],
            "sub_kI_B": sA["kI_B"],
        }
        comparisons.append(row)
        log(
            "ell/H=%.1f  gA_fem=%.16f  gA_sub=%.16f  |d|=%.3e"
            % (eoh, gA, sA["gamma_A"], dA_sub)
        )
        log(
            "         gB_fem=%.16f  gB_sub=%.16f  |d|=%.3e"
            % (gB, sA["gamma_B"], dB_sub)
        )
        log(
            "         kA_fem=%.16f  kA_sub=%.16f  |d|=%.3e"
            % (kA, sA["kI_A"], dKA_sub)
        )
        log(
            "         kB_fem=%.16f  kB_sub=%.16f  |d|=%.3e"
            % (kB, sA["kI_B"], dKB_sub)
        )
        log(
            "         int v' = %.12f  (NOT zero; ANSWER.md linear-term claim)"
            % lin
        )
        log("         natural-BC diagnostic vp0=%.6f vpL=%.6f (target -1)" % (vp0, vpL))

    log("")
    log("[FORMULA] reviewer tanh vs submission (not used inside FEM)")
    n_formula_agree = 0
    for eoh in SAMPLES:
        HoL = 1.0 / eoh
        gA_f = math.tanh(HoL / 2.0)
        gB_f = 2.0 * math.tanh(HoL / 2.0) - math.tanh(HoL)
        kA_f = math.sqrt(ONE_MINUS_NU2 * gA_f)
        kB_f = math.sqrt(ONE_MINUS_NU2 * gB_f)
        sA = sub_rows[eoh]
        for name, a, b in (
            ("gamma_A", gA_f, sA["gamma_A"]),
            ("gamma_B", gB_f, sA["gamma_B"]),
            ("kI_A", kA_f, sA["kI_A"]),
            ("kI_B", kB_f, sA["kI_B"]),
        ):
            d = abs(a - b)
            ok = d < 1e-12
            if ok:
                n_formula_agree += 1
            log("  ell/H=%.1f %s formula=%.16f sub=%.16f |d|=%.3e" % (eoh, name, a, b, d))
    log("formula vs submission agreement at 1e-12: %d/16" % n_formula_agree)
    log("FEM vs submission agreement at 1e-12: %d/16 (expect FEM truncation > 1e-12)" % n_agree_1e12)
    log("FEM vs submission agreement at 1e-6: %d/16" % n_agree_1e6)

    log("")
    log("[LIMITS] small-H series from FEM vs claimed expansions")
    for HoL in (0.02, 0.05):
        _, I_H, _, _, _ = solve_p1(HoL, 2000)
        _, I_2H, _, _, _ = solve_p1(2.0 * HoL, 2000)
        gA = -I_H
        gB = -2.0 * I_H + I_2H
        gA_s = HoL / 2.0
        gB_s = (HoL ** 3) / 4.0
        log(
            "  H/ell=%.3f  gA_fem=%.6e ~ H/2ell=%.6e  rel=%.3e"
            % (HoL, gA, gA_s, abs(gA - gA_s) / gA_s)
        )
        log(
            "             gB_fem=%.6e ~ H^3/4ell^3=%.6e  rel=%.3e"
            % (gB, gB_s, abs(gB - gB_s) / gB_s)
        )

    log("")
    log("[LIMITS] large-H FEM vs isolated gamma->1")
    for HoL in (8.0, 16.0):
        _, I_H, _, _, _ = solve_p1(HoL, 4000)
        _, I_2H, _, _, _ = solve_p1(2.0 * HoL, 4000)
        gA = -I_H
        gB = -2.0 * I_H + I_2H
        log("  H/ell=%.1f  gA_fem=%.10f  gB_fem=%.10f" % (HoL, gA, gB))

    log("")
    log("[CONTROL] C1 omit foundation energy at ell/H samples (should NOT match)")
    for eoh in SAMPLES:
        HoL = 1.0 / eoh
        g_wrong = control_omit_foundation(HoL, 2000)
        g_sub = sub_rows[eoh]["gamma_A"]
        log(
            "  ell/H=%.1f  C1_gamma_A=%.6f  sub=%.6f  |d|=%.3e"
            % (eoh, g_wrong, g_sub, abs(g_wrong - g_sub))
        )

    log("")
    log("[CONTROL] C2 incremental traction (Delta_sigma n=0) at L=1")
    g2, vmax = control_incremental_traction(1.0, 2000)
    log("  gamma=%.6e  max|v|=%.6e  (expect ~0; detects wrong face BC)" % (g2, vmax))

    log("")
    log("[LINEAR TERM] int h sigma0 eps_gg is NOT zero")
    log("  ANSWER.md L88-90 claims it integrates to zero; FEM int v' = v(L)-v(0)")
    log("  For odd opening about midplane, v(L)=-v(0) so int v' = -2 v(0) != 0")
    log("  compute.py still includes the term in W_dens; closed form remains tanh.")

    log("")
    log("[MODE] G=KI^2/E => kI=sqrt((1-nu^2) gamma); KII=0 by reflection.")
    log("[BOOKKEEPING] Case A: period H, 1 moving tip, divide by h.")
    log("[BOOKKEEPING] Case B: period 2H, 1 new moving tip; ahead S=2H, behind 2x S=H.")

    payload = {
        "n_el_production": n_el,
        "formula_vs_submission_1e12": n_formula_agree,
        "fem_vs_submission_1e12": n_agree_1e12,
        "fem_vs_submission_1e6": n_agree_1e6,
        "fem_conv_ratio_640_1280": ratio,
        "comparisons": comparisons,
        "linear_term_nonzero": True,
        "controls": {
            "C2_gamma": g2,
            "C2_max_v": vmax,
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(payload, f, indent=2)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines) + "\n")
    log("")
    log("wrote %s and %s" % (OUT_TXT, OUT_JSON))


if __name__ == "__main__":
    main()
