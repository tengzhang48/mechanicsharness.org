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

from crack_tip_checks import check_field, williams_exponent_symmetric

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


def williams_notch(alpha, n_theta=4001):
    p = williams_exponent_symmetric(alpha)
    lam = 1 - p
    s = lam + 1
    th = np.linspace(0, alpha, n_theta)
    B = -np.cos((lam + 1) * alpha) / np.cos((lam - 1) * alpha)
    F = np.cos((lam + 1) * th) + B * np.cos((lam - 1) * th)
    Fp = -(lam + 1) * np.sin((lam + 1) * th) - B * (lam - 1) * np.sin((lam - 1) * th)
    Fpp = -(lam + 1)**2 * np.cos((lam + 1) * th) - B * (lam - 1)**2 * np.cos((lam - 1) * th)
    return th, Fpp + s * F, s * (s - 1) * F, -(s - 1) * Fp, p


def hrr_field(n, plane_strain):
    sol = hrr.solve_angular(n, plane_strain=plane_strain)
    return sol["theta"], sol["sr"], sol["st"], sol["srt"]


def power_notch(n=4, alpha=2 * np.pi / 3):
    spec = importlib.util.spec_from_file_location("t2ref", T2_REF)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    sol, info = m.shooting(n, alpha, (1.8438247397, -1.2630789963))
    th = np.linspace(0, alpha, 3001)
    f, fp, fpp, _ = sol.sol(th)
    s = info["s"]
    return th, fpp + s * f, s * (s - 1) * f, -(s - 1) * fp, info["p"]


def test_williams_crack_passes_and_J_is_K2_over_E_prime():
    rep = check_field(*williams_crack(), p=0.5, law="linear", nu=0.3, plane="strain")
    assert failed(rep) == set()
    assert abs(by_name(rep)["J of the field (report)"]["value"] - (1 - 0.3**2)) < 1e-6


def test_williams_notch_passes():
    th, a, b, c, p = williams_notch(2 * np.pi / 3)
    rep = check_field(th, a, b, c, p=p, law="linear", alpha=2 * np.pi / 3)
    assert failed(rep) == set()
    assert by_name(rep)["linear notch exponent (Williams, mode I)"]["passed"] is True


def test_williams_notch_exponent_known_value():
    # 120 degree notch (faces at +-120 deg): lambda = 0.6157 (also the probe builder's linear check)
    assert abs((1 - williams_exponent_symmetric(2 * np.pi / 3)) - 0.6157) < 1e-4


@pytest.mark.parametrize("n,pe,I_table", [(3, False, 3.86), (3, True, 5.51), (9, True, 4.60)])
def test_hrr_fields_pass_and_J_equals_table_I(n, pe, I_table):
    rep = check_field(*hrr_field(n, pe), p=1 / (n + 1), n=n, plane="strain" if pe else "stress")
    assert failed(rep) == set()
    assert abs(by_name(rep)["J of the field (report)"]["value"] - I_table) / I_table < 0.005


def test_power_law_notch_passes():
    th, a, b, c, p = power_notch()
    rep = check_field(th, a, b, c, p=p, n=4, alpha=2 * np.pi / 3, plane="strain")
    assert failed(rep) == set()
    assert by_name(rep)["crack exponent fixed by J"]["applies"] is False


# ------------------------------------------------------------------ planted errors
def test_wrong_exponent_fails_equilibrium_and_crack_exponent():
    rep = check_field(*hrr_field(3, True), p=0.26, n=3, plane="strain")
    assert {"equilibrium (angular ODEs, declared p)", "crack exponent fixed by J"} <= failed(rep)


def test_notch_with_inconsistent_exponent_fails_equilibrium():
    # the HRR-notch probe failure (grader R23): equilibrium used with the wrong radial exponent
    th, a, b, c, p = williams_notch(2 * np.pi / 3)
    rep = check_field(th, a, b, c, p=p + 0.02, law="linear", alpha=2 * np.pi / 3)
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
    rep = check_field(th[keep], a[keep], b[keep], c[keep], p=0.5, law="linear", alpha=th[keep][-1])
    assert "faces traction-free" in failed(rep)


def test_broken_symmetry_is_flagged():
    th, a, b, c = williams_crack()
    rep = check_field(th, a + 0.01 * np.max(np.abs(a)) * np.sin(th), b, c, p=0.5, law="linear")
    assert "mode I symmetry" in failed(rep)


def test_infinite_energy_exponent_is_flagged():
    th, a, b, c = hrr_field(3, True)
    rep = check_field(th, a, b, c, p=0.6, n=3, plane="strain")
    assert "finite energy near the tip" in failed(rep)
