#!/usr/bin/env python3
"""Independent deciding-quantity checks for the Xia-Hutchinson strip energies.

Does NOT import or call the author's analytic functions. Solves a variational
discretization of the PROBLEM.md full-width strip functional whose natural BCs
are the total-traction face conditions.

Single command: python3 check_independent.py
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NU = 0.3
ONE_MINUS_NU2 = 1.0 - NU * NU
SAMPLES = [0.5, 1.0, 1.5, 2.0]
OUT_TXT = os.path.join(HERE, "check_independent_output.txt")
REF_JSON = os.path.join(HERE, "reference_submission", "reference.json")


def log(lines, s=""):
    print(s)
    lines.append(s)


# ---------------------------------------------------------------------------
# Independent closed forms (derived in DERIVATION.md; not copied from author
# function names). Used only as a formula check, not inside the FEM.
# ---------------------------------------------------------------------------

def gamma_A_formula(H_over_ell):
    return math.tanh(H_over_ell / 2.0)


def gamma_B_formula(H_over_ell):
    return 2.0 * math.tanh(H_over_ell / 2.0) - math.tanh(H_over_ell)


def kI_formula(gamma):
    return math.sqrt(ONE_MINUS_NU2 * gamma)


# ---------------------------------------------------------------------------
# Variational P1 FEM on the FULL period (0, L), L = S/ell.
# Discrete energy:
#   I[v] = v(L) - v(0) + int_0^L 1/2 (v'^2 + v^2) d eta
# Stationarity => v'' = v in the interior and natural BCs v'(0)=v'(L)=-1
# (total traction sigma22=0 on both crack faces). No essential slope BCs.
# ---------------------------------------------------------------------------

def assemble_p1(L, n_el):
    n = n_el + 1
    h = L / n_el
    K = np.zeros((n, n))
    ke_diff = (1.0 / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])
    ke_mass = (h / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
    ke = ke_diff + ke_mass
    for e in range(n_el):
        K[e : e + 2, e : e + 2] += ke
    rhs = np.zeros(n)
    rhs[0] = 1.0
    rhs[-1] = -1.0
    return K, rhs, h


def solve_fullwidth_p1(L, n_el):
    K, rhs, h = assemble_p1(L, n_el)
    v = np.linalg.solve(K, rhs)
    I = (v[-1] - v[0]) + 0.5 * float(v @ (K @ v))
    eta = np.linspace(0.0, L, n_el + 1)
    return eta, v, I, h


def gamma_from_I(I):
    return -I


def face_slopes_p1(v, h):
    # one-sided first differences as a diagnostic of the natural BC
    vp0 = (v[1] - v[0]) / h
    vpL = (v[-1] - v[-2]) / h
    return vp0, vpL


# ---------------------------------------------------------------------------
# Independent second-order variational FD (ghost-free energy on nodes).
# I_h = v_N - v_0 + sum_i 1/2 ( (v_{i+1}-v_i)/h )^2 h + trapezoid 1/2 v^2
# This is the same P1 energy, so it is a cross-check of the assembler.
# A genuinely different route: collocation-free weak FD with trapezoid mass
# and natural (Neumann) rows from energy stationarity, using a lumped mass.
# ---------------------------------------------------------------------------

def solve_fullwidth_lumped(L, n_el):
    """Energy with lumped (trapezoid/vertex) mass: different from consistent P1.

    I = v_N - v_0 + (h/2) sum_i ((v_{i+1}-v_i)/h)^2 * 2/2 wait:
    int v'^2/2 ~ sum_e (dv)^2 /(2h)
    int v^2/2 ~ (h/2) * (v0^2/2 + vN^2/2 + sum_{interior} v_i^2)  [trapezoid]
    """
    n = n_el + 1
    h = L / n_el
    # stiffness from int v'^2  (so energy uses 1/2 v^T A v with A = int (v' v' + v v))
    A = np.zeros((n, n))
    for e in range(n_el):
        A[e, e] += 1.0 / h
        A[e, e + 1] -= 1.0 / h
        A[e + 1, e] -= 1.0 / h
        A[e + 1, e + 1] += 1.0 / h
    # lumped mass int v^2 : trapezoid weights
    w = np.full(n, h)
    w[0] = 0.5 * h
    w[-1] = 0.5 * h
    A += np.diag(w)
    rhs = np.zeros(n)
    rhs[0] = 1.0
    rhs[-1] = -1.0
    v = np.linalg.solve(A, rhs)
    I = (v[-1] - v[0]) + 0.5 * float(v @ (A @ v))
    return v, I


# ---------------------------------------------------------------------------
# Small-H series (independent Taylor), isolated limit, controls
# ---------------------------------------------------------------------------

def series_gamma_A(H_over_ell):
    x = H_over_ell / 2.0
    # tanh(x) = x - x^3/3 + 2 x^5/15
    return x - x**3 / 3.0 + 2.0 * x**5 / 15.0


def series_gamma_B(H_over_ell):
    # 2 tanh(x) - tanh(2x) with x = H/(2 ell) = 2 x^3 - 4 x^5 + ...
    x = H_over_ell / 2.0
    return 2.0 * x**3 - 4.0 * x**5


def control_C1_omit_foundation_exact(a):
    # I = int v' + (v')^2/2  using exact field (diagnostic of the error class,
    # not used as a production gamma). wrong gamma_A = -I
    # derived: 1.5 tanh(a) - (a/2) sech^2(a)
    sech2 = 1.0 / math.cosh(a) ** 2
    return 1.5 * math.tanh(a) - 0.5 * a * sech2


def control_C1_omit_foundation_fem(L, n_el):
    """Same FEM but drop the v^2 (foundation) term from the energy."""
    n = n_el + 1
    h = L / n_el
    K = np.zeros((n, n))
    ke_diff = (1.0 / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])
    for e in range(n_el):
        K[e : e + 2, e : e + 2] += ke_diff
    rhs = np.zeros(n)
    rhs[0] = 1.0
    rhs[-1] = -1.0
    # Without foundation, K is singular (constants). Use the equilibrium field
    # from the FULL energy, then evaluate the truncated functional (author C1).
    Kfull, rhsfull, _ = assemble_p1(L, n_el)
    v = np.linalg.solve(Kfull, rhsfull)
    I = (v[-1] - v[0]) + 0.5 * float(v @ (K @ v))  # K is diffusion only
    return -I


def control_C2_incremental_fem(L, n_el):
    """Natural BC v'=0 (incremental traction) instead of v'=-1: drop linear term."""
    K, _, _ = assemble_p1(L, n_el)
    rhs = np.zeros(n_el + 1)
    # K v = 0 => v = 0 (K SPD)
    v = np.linalg.solve(K, rhs)
    I = (v[-1] - v[0]) + 0.5 * float(v @ (K @ v))
    return -I, float(np.max(np.abs(v)))


def main():
    lines = []
    log(lines, "=" * 88)
    log(lines, "Independent checker: full-width variational strip (PROBLEM.md energy)")
    log(lines, "nu=%.1f, 1-nu^2=%.2f, ell=sqrt(Ebar h/k)" % (NU, ONE_MINUS_NU2))
    log(lines, "Natural BCs: v'(0)=v'(L)=-1 (total traction). No author functions.")
    log(lines, "=" * 88)

    with open(REF_JSON) as f:
        ref = json.load(f)
    ref_rows = {r["ell_over_H"]: r for r in ref["rows"]}

    # --- FEM convergence at one spacing (S/ell = 1, i.e. ell/H=1 Case A) ---
    log(lines, "")
    log(lines, "[FEM-conv] full-width P1 energy at L=S/ell=1 (a=0.5), vs tanh(0.5)")
    target = math.tanh(0.5)
    conv = []
    for n_el in (40, 80, 160, 320, 640, 1280, 2560):
        _, v, I, h = solve_fullwidth_p1(1.0, n_el)
        g = gamma_from_I(I)
        err = abs(g - target)
        vp0, vpL = face_slopes_p1(v, h)
        conv.append((n_el, err, vp0, vpL, g))
        log(
            lines,
            "  n_el=%4d  gamma=%.12f  |err|=%.3e  vp(0)=%.6f vp(L)=%.6f"
            % (n_el, g, err, vp0, vpL),
        )
    # P1 energy for this 1D Helmholtz problem is O(h^2)
    e640 = conv[-3][1]
    e1280 = conv[-2][1]
    ratio = e640 / e1280 if e1280 else float("inf")
    log(lines, "  observed err(640)/err(1280)=%.3f (expect ~4 for O(h^2))" % ratio)

    # --- Production FEM table vs formula vs reference.json ---
    log(lines, "")
    log(lines, "[MAIN] production full-width P1 (n_el=4000) and lumped-mass (n_el=8000)")
    log(
        lines,
        "%8s %8s %16s %16s %16s %16s %12s %12s"
        % ("ell/H", "H/ell", "gA_fem", "gB_fem", "gA_ref", "gB_ref", "dA", "dB"),
    )
    n_el = 4000
    n_el_lump = 8000
    comparisons = []
    for eoh in SAMPLES:
        HoL = 1.0 / eoh
        _, _, I_H, _ = solve_fullwidth_p1(HoL, n_el)
        _, _, I_2H, _ = solve_fullwidth_p1(2.0 * HoL, n_el)
        gA_fem = gamma_from_I(I_H)
        # Case B: one new tip per period 2H: -I_behind_over_2H + I_ahead
        # behind = two H-cells = 2 I(H); ahead = I(2H)
        # G_B * (Ebar/(sigma0^2 ell)) = - (2 I(H) - I(2H)) wait:
        # I(S) already = int over that period. Pi ~ I.
        # G_B h / (ell h sigma0^2/Ebar) = I(2H) wait no:
        # energy release = Pi_ahead - Pi_behind, Pi ~ I (same prefactor).
        # Pi_ahead over width 2H corresponds to I(2H)
        # Pi_behind over width 2H corresponds to 2*I(H)
        # gamma_B = I(2H) - 2*I(H)   but I is negative, gamma_B = -2 I(H) + I(2H)
        # since gamma_A = -I(H), gamma_B = 2 gA(H) - gA(2H) = -2 I(H) + I(2H)
        gB_fem = -2.0 * I_H + I_2H
        _, I_H_l = solve_fullwidth_lumped(HoL, n_el_lump)
        _, I_2H_l = solve_fullwidth_lumped(2.0 * HoL, n_el_lump)
        gA_lump = -I_H_l
        gB_lump = -2.0 * I_H_l + I_2H_l
        gA_f = gamma_A_formula(HoL)
        gB_f = gamma_B_formula(HoL)
        kA = kI_formula(gA_f)
        kB = kI_formula(gB_f)
        r = ref_rows[eoh]
        dA = gA_fem - r["gamma_A"]
        dB = gB_fem - r["gamma_B"]
        dA_f = gA_f - r["gamma_A"]
        dB_f = gB_f - r["gamma_B"]
        dkA = kA - r["kI_A"]
        dkB = kB - r["kI_B"]
        log(
            lines,
            "%8.4f %8.6f %16.12f %16.12f %16.12f %16.12f %12.3e %12.3e"
            % (eoh, HoL, gA_fem, gB_fem, r["gamma_A"], r["gamma_B"], dA, dB),
        )
        comparisons.append(
            dict(
                ell_over_H=eoh,
                H_over_ell=HoL,
                gamma_A_formula=gA_f,
                gamma_B_formula=gB_f,
                kI_A_formula=kA,
                kI_B_formula=kB,
                gamma_A_fem_p1=gA_fem,
                gamma_B_fem_p1=gB_fem,
                gamma_A_fem_lumped=gA_lump,
                gamma_B_fem_lumped=gB_lump,
                gamma_A_ref=r["gamma_A"],
                gamma_B_ref=r["gamma_B"],
                kI_A_ref=r["kI_A"],
                kI_B_ref=r["kI_B"],
                d_gamma_A_fem_vs_ref=dA,
                d_gamma_B_fem_vs_ref=dB,
                d_gamma_A_formula_vs_ref=dA_f,
                d_gamma_B_formula_vs_ref=dB_f,
                d_kI_A=dkA,
                d_kI_B=dkB,
            )
        )

    max_fem = max(abs(c["d_gamma_A_fem_vs_ref"]) for c in comparisons + [])
    max_fem = max(max_fem, max(abs(c["d_gamma_B_fem_vs_ref"]) for c in comparisons))
    max_form = max(abs(c["d_gamma_A_formula_vs_ref"]) for c in comparisons)
    max_form = max(max_form, max(abs(c["d_gamma_B_formula_vs_ref"]) for c in comparisons))
    max_k = max(abs(c["d_kI_A"]) for c in comparisons)
    max_k = max(max_k, max(abs(c["d_kI_B"]) for c in comparisons))
    max_lump = max(
        max(abs(c["gamma_A_fem_lumped"] - c["gamma_A_ref"]) for c in comparisons),
        max(abs(c["gamma_B_fem_lumped"] - c["gamma_B_ref"]) for c in comparisons),
    )
    log(lines, "")
    log(lines, "[MAIN] max |P1 FEM - ref gamma| = %.3e" % max_fem)
    log(lines, "[MAIN] max |lumped FEM - ref gamma| = %.3e" % max_lump)
    log(lines, "[MAIN] max |formula - ref gamma| = %.3e" % max_form)
    log(lines, "[MAIN] max |kI formula - ref| = %.3e" % max_k)

    # --- natural-BC diagnostic ---
    log(lines, "")
    log(lines, "[BC] P1 face slopes should approach -1 (total traction) as h->0")
    _, v, _, h = solve_fullwidth_p1(1.0, 2560)
    vp0, vpL = face_slopes_p1(v, h)
    log(lines, "  L=1, n_el=2560: vp(0)=%.8f  vp(L)=%.8f  (target -1)" % (vp0, vpL))

    # --- limits ---
    log(lines, "")
    log(lines, "[LIMITS] small-H series vs formula; isolated H/ell=20 vs 1")
    for HoL in (0.05, 0.1, 0.2):
        gA = gamma_A_formula(HoL)
        gB = gamma_B_formula(HoL)
        sA = series_gamma_A(HoL)
        sB = series_gamma_B(HoL)
        log(
            lines,
            "  H/ell=%.3f  gA=%.8e series=%.8e d=%.2e ; gB=%.8e series=%.8e d=%.2e"
            % (HoL, gA, sA, gA - sA, gB, sB, gB - sB),
        )
        # sigma_c scalings from leading terms
        # gamma_A ~ H/(2 ell); gamma_B ~ (H/ell)^3/4
        leadA = HoL / 2.0
        leadB = (HoL**3) / 4.0
        log(
            lines,
            "           lead gA=H/(2ell)=%.8e rel=%.3e ; lead gB=(H/ell)^3/4=%.8e rel=%.3e"
            % (leadA, (gA - leadA) / leadA, leadB, (gB - leadB) / leadB),
        )
    gA_inf = gamma_A_formula(20.0)
    gB_inf = gamma_B_formula(20.0)
    kI_iso = math.sqrt(ONE_MINUS_NU2)
    log(
        lines,
        "  H/ell=20: gA=%.10f gB=%.10f (target 1); kI(gamma=1)=%.10f (Eq27 sqrt(1-nu^2))"
        % (gA_inf, gB_inf, kI_iso),
    )

    # --- controls (independent) ---
    log(lines, "")
    log(lines, "[CONTROLS] independent reconstructions of author's wrong variants")
    c1_fem = []
    c1_exact = []
    for eoh in SAMPLES:
        HoL = 1.0 / eoh
        a = HoL / 2.0
        c1_exact.append(control_C1_omit_foundation_exact(a))
        c1_fem.append(control_C1_omit_foundation_fem(HoL, 2000))
    log(lines, "  C1 omit foundation (exact-field truncated energy) gamma_A_wrong=%s" %
        ["%.6f" % x for x in c1_exact])
    log(lines, "  C1 omit foundation (FEM truncated energy)          gamma_A_wrong=%s" %
        ["%.6f" % x for x in c1_fem])
    log(lines, "  C1 vs author saved %s" %
        ["%.6f" % x for x in ref["meta"]["controls_deliberately_wrong"]["C1_omit_foundation_energy_gamma_A_wrong"]])

    gC2, maxv = control_C2_incremental_fem(1.0, 400)
    log(lines, "  C2 incremental traction FEM: gamma_A=%.3e max|v|=%.3e (expect 0)" % (gC2, maxv))

    c3 = []
    for eoh in SAMPLES:
        HoL = 1.0 / eoh
        # swapped: 2 G'(2H) - G'(H)
        c3.append(2.0 * gamma_A_formula(2.0 * HoL) - gamma_A_formula(HoL))
    log(lines, "  C3 swapped spacings gamma_B_wrong=%s" % ["%.6f" % x for x in c3])
    log(lines, "  C3 vs author saved %s" %
        ["%.6f" % x for x in ref["meta"]["controls_deliberately_wrong"]["C3_swapped_spacings_gamma_B_wrong"]])

    # --- mode relation sanity ---
    log(lines, "")
    log(lines, "[MODE] G=KI^2/E, Ebar=E/(1-nu^2) => kI=sqrt((1-nu^2) gamma); KII=0 by symmetry")
    log(lines, "  sqrt(1-nu^2)=%.12f  (isolated-crack Eq. 27 target)" % kI_iso)

    # --- Eq 11-12 calibration identities (algebraic) ---
    log(lines, "")
    log(lines, "[CAL] isolated G_model = sigma0^2 ell/Ebar ; Beuth G = (pi/2)(h sigma0^2/Ebar) g")
    log(lines, "      matching => ell = (pi/2) g(alpha,beta) h   (XH2000 Eq. 12)")
    log(lines, "      no mismatch: paper states ell=1.98 h (not re-derived from Beuth tables here)")

    summary = dict(
        max_abs_d_gamma_fem_p1_vs_ref=max_fem,
        max_abs_d_gamma_fem_lumped_vs_ref=max_lump,
        max_abs_d_gamma_formula_vs_ref=max_form,
        max_abs_d_kI_formula_vs_ref=max_k,
        fem_conv_err_n40_to_n2560=[c[1] for c in conv],
        fem_conv_ratio_640_over_1280=ratio,
        face_slopes_n2560_L1={"vp0": vp0, "vpL": vpL},
        isolated_kI_Eq27=kI_iso,
        H20_gamma_A=gA_inf,
        H20_gamma_B=gB_inf,
        C1_exact=c1_exact,
        C1_fem=c1_fem,
        C2_gamma=gC2,
        C2_maxv=maxv,
        C3=c3,
        rows=comparisons,
        shared_assumptions=[
            "plane-stress Winkler strip; u1=0 (mode I) so the check is 1D wake energy, not a 2D tip field",
            "PROBLEM.md incremental functional; total (not incremental) face traction",
            "steady translating-tip conversion of far-ahead strip to far-behind wake",
            "Case A: one tip per period H; Case B: one new tip per period 2H",
            "G=(KI^2+KII^2)/E with KII=0",
        ],
        what_this_detects=[
            "wrong energy density (missing foundation or prestress term)",
            "incremental vs total traction on faces",
            "wrong period / tip count / ahead-behind assignment",
            "E vs Ebar mix-up in kI",
            "does NOT detect 2D tip-field or wavenumber-dependent substrate errors",
        ],
    )
    out_json = os.path.join(HERE, "check_independent_results.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)
    log(lines, "")
    log(lines, "[OUT] wrote %s" % out_json)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines) + "\n")
    log(lines, "[OUT] wrote %s" % OUT_TXT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
