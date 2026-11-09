"""Validation of crack_tip_checks: green on correct fields, red on planted errors (log P240).

Correct fields come from closed forms (Williams crack and notch), from the frozen HRR v1 code
(knowledge_v1/P015_hutchinson_1968) and from the HRR-notch probe's builder reference (Codex,
P207). The checks themselves import none of these.
Run: python3 -m pytest -q test_checks.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from crack_tip_checks import check_field, williams_exponent

PILOTS = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PILOTS / "corpus/knowledge_v1/P015_hutchinson_1968"))
import hrr  # noqa: E402

T2_REF = Path("<probe-runs>/P015-build-codex/T2/reference.py")


def by_name(report):
    return {r["check"]: r for r in report}


def failed(report):
    return {r["check"] for r in report if r["passed"] is False}


# ------------------------------------------------------------------ correct fields
def williams_crack(n_theta=4001):
    th = np.linspace(-np.pi, np.pi, n_theta)
    k = 1 / np.sqrt(2 * np.pi)
    return (th, k * (1.25 * np.cos(th / 2) - 0.25 * np.cos(1.5 * th)),
            k * (0.75 * np.cos(th / 2) + 0.25 * np.cos(1.5 * th)),
            k * (0.25 * np.sin(th / 2) + 0.25 * np.sin(1.5 * th)))


def williams_notch(beta, n_theta=4001):
    p = williams_exponent(beta)
    lam = 1 - p
    s = lam + 1
    th = np.linspace(0, beta, n_theta)
    B = -np.cos((lam + 1) * beta) / np.cos((lam - 1) * beta)
    F = np.cos((lam + 1) * th) + B * np.cos((lam - 1) * th)
    Fp = -(lam + 1) * np.sin((lam + 1) * th) - B * (lam - 1) * np.sin((lam - 1) * th)
    Fpp = -(lam + 1)**2 * np.cos((lam + 1) * th) - B * (lam - 1)**2 * np.cos((lam - 1) * th)
    return th, Fpp + s * F, s * (s - 1) * F, -(s - 1) * Fp, p


def hrr_field(n, plane_strain):
    sol = hrr.solve_angular(n, plane_strain=plane_strain)
    return sol["theta"], sol["sr"], sol["st"], sol["srt"]


def power_notch(n=4, beta=2 * np.pi / 3):
    spec = importlib.util.spec_from_file_location("t2ref", T2_REF)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    sol, info = m.shooting(n, beta, (1.8438247397, -1.2630789963))
    th = np.linspace(0, beta, 3001)
    f, fp, fpp, _ = sol.sol(th)
    s = info["s"]
    return th, fpp + s * f, s * (s - 1) * f, -(s - 1) * fp, info["p"]


def test_williams_crack_passes_and_J_is_K2_over_E_prime():
    rep = check_field(*williams_crack(), p=0.5, law="linear", nu=0.3, plane="strain")
    assert failed(rep) == set()
    assert abs(by_name(rep)["J of the field (report)"]["value"] - (1 - 0.3**2)) < 1e-6


def test_williams_notch_passes():
    th, a, b, c, p = williams_notch(2 * np.pi / 3)
    rep = check_field(th, a, b, c, p=p, law="linear", beta=2 * np.pi / 3)
    assert failed(rep) == set()
    assert by_name(rep)["linear notch exponent (Williams, mode I)"]["passed"] is True


def test_williams_notch_exponent_known_value():
    # 120 degree notch (faces at +-120 deg): lambda = 0.6157 (also the probe builder's linear check)
    assert abs((1 - williams_exponent(2 * np.pi / 3)) - 0.6157) < 1e-4


@pytest.mark.parametrize("n,pe,I_table", [(3, False, 3.86), (3, True, 5.51), (9, True, 4.60)])
def test_hrr_fields_pass_and_J_equals_table_I(n, pe, I_table):
    rep = check_field(*hrr_field(n, pe), p=1 / (n + 1), n=n, plane="strain" if pe else "stress")
    assert failed(rep) == set()
    assert abs(by_name(rep)["J of the field (report)"]["value"] - I_table) / I_table < 0.005


def test_power_law_notch_passes():
    th, a, b, c, p = power_notch()
    rep = check_field(th, a, b, c, p=p, n=4, beta=2 * np.pi / 3, plane="strain")
    assert failed(rep) == set()
    assert by_name(rep)["crack exponent fixed by J"]["applies"] is False


# ------------------------------------------------------------------ planted errors
def test_wrong_exponent_fails_equilibrium_and_crack_exponent():
    rep = check_field(*hrr_field(3, True), p=0.26, n=3, plane="strain")
    assert {"equilibrium (angular ODEs, declared p)", "crack exponent fixed by J"} <= failed(rep)


def test_notch_with_inconsistent_exponent_fails_equilibrium():
    # the HRR-notch probe failure (grader R23): equilibrium used with the wrong radial exponent
    th, a, b, c, p = williams_notch(2 * np.pi / 3)
    rep = check_field(th, a, b, c, p=p + 0.02, law="linear", beta=2 * np.pi / 3)
    assert "equilibrium (angular ODEs, declared p)" in failed(rep)
    assert "equilibrium (Cartesian, finite differences)" in failed(rep)


def test_scaled_shear_component_fails_equilibrium():
    th, a, b, c = hrr_field(3, True)
    rep = check_field(th, a, b, 1.05 * c, p=0.25, n=3, plane="strain")
    assert "equilibrium (angular ODEs, declared p)" in failed(rep)


def test_wrong_hardening_exponent_fails_compatibility():
    rep = check_field(*hrr_field(3, True), p=0.25, n=5, plane="strain")
    assert "compatibility (displacements from strains, shear relation)" in failed(rep)


def test_plane_stress_field_declared_plane_strain_fails_compatibility():
    rep = check_field(*hrr_field(3, False), p=0.25, n=3, plane="strain")
    assert "compatibility (displacements from strains, shear relation)" in failed(rep)


def test_reported_I_from_dropped_basis_derivatives_is_flagged():
    # the HRR unit's original error: I = 5.65 instead of 3.86 (plane stress, n = 3)
    rep = check_field(*hrr_field(3, False), p=0.25, n=3, plane="stress", reported_J=5.6476)
    assert "reported J (or I) matches the field's own J" in failed(rep)
    rep_ok = check_field(*hrr_field(3, False), p=0.25, n=3, plane="stress", reported_J=3.855)
    assert "reported J (or I) matches the field's own J" not in failed(rep_ok)


def test_faces_not_traction_free_are_flagged():
    th, a, b, c = williams_crack()
    keep = np.abs(th) <= 0.9 * np.pi + 1e-12
    rep = check_field(th[keep], a[keep], b[keep], c[keep], p=0.5, law="linear", beta=th[keep][-1])
    assert "faces traction-free" in failed(rep)


def test_broken_symmetry_is_flagged():
    th, a, b, c = williams_crack()
    rep = check_field(th, a + 0.01 * np.max(np.abs(a)) * np.sin(th), b, c, p=0.5, law="linear")
    assert "mode I symmetry" in failed(rep)


def test_infinite_energy_exponent_is_flagged():
    th, a, b, c = hrr_field(3, True)
    rep = check_field(th, a, b, c, p=0.6, n=3, plane="strain")
    assert "finite energy near the tip" in failed(rep)


# ------------------------------------------------------------------ cases from Grok's review (P244)
def test_t_stress_passes_as_a_non_leading_term():
    # sigma_xx = T, p = 0: correct, traction-free, mode I; previously failed the relative Cartesian checks
    th = np.linspace(-np.pi, np.pi, 4001)
    rep = check_field(th, np.cos(th)**2, np.sin(th)**2, -np.sin(th) * np.cos(th), p=0.0, law="linear",
                      leading=False)
    assert failed(rep) == set()


def test_williams_mode_II_crack_passes():
    th = np.linspace(-np.pi, np.pi, 4001)
    k = 1 / np.sqrt(2 * np.pi)
    srr = k * (-1.25 * np.sin(th / 2) + 0.75 * np.sin(1.5 * th))
    stt = k * (-0.75 * np.sin(th / 2) - 0.75 * np.sin(1.5 * th))
    srt = k * (0.25 * np.cos(th / 2) + 0.75 * np.cos(1.5 * th))
    rep = check_field(th, srr, stt, srt, p=0.5, law="linear", mode="II")
    assert failed(rep) == set()


def test_90_degree_wedge_is_non_singular_and_does_not_crash():
    # beta = pi/2 (a half-plane edge): lambda = 1, p = 0; the field is a uniform stress along the edge
    assert abs(williams_exponent(np.pi / 2)) < 1e-9
    th = np.linspace(-np.pi / 2, np.pi / 2, 4001)
    rep = check_field(th, 4 * np.sin(th)**2, 4 * np.cos(th)**2, 4 * np.sin(th) * np.cos(th), p=0.0,
                      law="linear", beta=np.pi / 2)
    assert failed(rep) == set()


def test_logarithmic_displacement_case_does_not_crash():
    th, a, b, c = williams_crack()
    rep = check_field(th, a, b, c, p=1.0, law="linear")
    assert by_name(rep)["displacement and J checks"]["applies"] is False
    assert "finite energy near the tip" in failed(rep)


def test_power_law_notch_exponent_error_of_1e4_is_caught():
    th, a, b, c, p = power_notch()
    rep = check_field(th, a, b, c, p=p + 1e-4, n=4, beta=2 * np.pi / 3, plane="strain")
    assert "equilibrium (angular ODEs, declared p)" in failed(rep)


def test_williams_mode_II_notch_150_degrees():
    beta = 5 * np.pi / 6
    p = williams_exponent(beta, "II")
    assert abs(p - 0.269) < 1e-3          # Grok's independent value
    lam = 1 - p
    s = lam + 1
    th = np.linspace(0, beta, 4001)
    B = -np.sin((lam + 1) * beta) / np.sin((lam - 1) * beta)
    F = np.sin((lam + 1) * th) + B * np.sin((lam - 1) * th)
    Fp = (lam + 1) * np.cos((lam + 1) * th) + B * (lam - 1) * np.cos((lam - 1) * th)
    Fpp = -(lam + 1)**2 * np.sin((lam + 1) * th) - B * (lam - 1)**2 * np.sin((lam - 1) * th)
    rep = check_field(th, Fpp + s * F, s * (s - 1) * F, -(s - 1) * Fp, p=p, law="linear", beta=beta, mode="II")
    assert failed(rep) == set()
    assert by_name(rep)["linear notch exponent (Williams, mode II)"]["passed"] is True


def test_closing_crack_is_flagged_as_interpenetration():
    # negative K_I passed every check before (Grok): traction-free faces need the faces to open
    th, a, b, c = williams_crack()
    rep = check_field(th, -a, -b, -c, p=0.5, law="linear")
    assert failed(rep) == {"crack faces open (no interpenetration)"}
    assert "crack faces open (no interpenetration)" not in failed(check_field(th, a, b, c, p=0.5, law="linear"))


def test_crack_face_displacement_matches_textbook_value():
    # upper-face u_y = (K/2mu) sqrt(r/2pi) (kappa+1), plane strain, K = E = 1, nu = 0.3, r = 1
    import crack_tip_checks as C
    th, a, b, c = williams_crack()
    t, S = C._full_range(th, a, b, c, "I")
    Ur, Ut, _ = C._Field(t, S, 0.5, 1, "linear", "strain", 0.3).displacement_splines()
    mu, kappa = 1 / (2 * 1.3), 3 - 4 * 0.3
    assert abs(-Ut(np.pi) - (1 / (2 * mu)) * np.sqrt(1 / (2 * np.pi)) * (kappa + 1)) < 1e-6


def test_singularity_thresholds():
    # mode I singular only for a re-entrant corner (total angle > 180 deg); mode II only beyond the
    # angle where lambda = 1 is a double root of the mode II equation, tan(2 beta) = 2 beta (257.45 deg)
    from scipy.optimize import brentq
    crit = brentq(lambda x: np.tan(x) - x, 4.0, 4.6) / 2
    assert williams_exponent(np.radians(179) / 2, "I") < 0 < williams_exponent(np.radians(181) / 2, "I")
    assert williams_exponent(crit - 0.005, "II") < 0 < williams_exponent(crit + 0.005, "II")
