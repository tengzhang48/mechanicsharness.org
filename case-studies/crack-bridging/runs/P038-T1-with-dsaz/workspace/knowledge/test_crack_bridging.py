"""Tests for the reproduction of Bao & Suo (1992) crack-bridging results.

Run with:  python3 -m pytest

Tolerances are justified by measured convergence errors (see REPORT.md):
  * zone solver (semi-infinite kernel) is exact to ~1e-4 on Dugdale checks;
  * spectral solver (finite center crack) has a traction-jump singularity at
    the trailing edge of the bridging zone: measured errors ~0.3-0.5 % on
    opening, ~0.05 % on K at N = 600-1500.
"""

import numpy as np
import pytest

import crack_bridging as cb
from zone import ssb_zone_solve, ssb_coefficients, ssb_rcurve_c1


# ---------------------------------------------------------------------------
# Table 1 (p. 356)
# ---------------------------------------------------------------------------

def test_table1():
    t1 = cb.table1()
    # transcription verified against pages/page-02.png (300 dpi)
    assert t1['atomic bond']['sigma0'] == pytest.approx(1e10)
    assert t1['atomic bond']['delta0'] == pytest.approx(1e-10)
    # products as printed in the paper; the table is internally consistent
    assert t1['atomic bond']['sigma0_delta0_printed'] == pytest.approx(1.0)
    assert t1['fiber pull-out']['sigma0_delta0_printed'] == pytest.approx(1e4)
    assert t1['fiber cross-over']['sigma0_delta0_printed'] == pytest.approx(1e3)
    assert t1['atomic bond']['delta0_E_over_sigma0_printed'] == pytest.approx(1e-9)
    assert t1['fiber pull-out']['delta0_E_over_sigma0_printed'] == pytest.approx(1e-3)
    assert t1['fiber cross-over']['delta0_E_over_sigma0_printed'] == pytest.approx(1.0)
    # printed products and lengths follow exactly from the printed sigma0, delta0
    for name in t1:
        assert t1[name]['sigma0_delta0_computed'] == \
            pytest.approx(t1[name]['sigma0_delta0_printed'], rel=1e-12)
        assert t1[name]['delta0_E_over_sigma0'] == \
            pytest.approx(t1[name]['delta0_E_over_sigma0_printed'], rel=1e-12)


# ---------------------------------------------------------------------------
# Kernel validation (closed forms, machine precision)
# ---------------------------------------------------------------------------

def test_kernel_remote_only():
    sol = cb.solve_bridged(1e-12, s=0.0, sigma0=0.0, remote=1.0, N=100)
    assert sol['K_tip'] == pytest.approx(np.sqrt(np.pi), rel=1e-12)
    assert sol['delta_at'](np.array([np.pi / 2]))[0] == pytest.approx(4.0, rel=1e-12)


def test_kernel_whole_crack_closure():
    sol = cb.solve_bridged(1.0, s=0.0, sigma0=1.0, remote=0.0, N=100)
    assert sol['K_tip'] == pytest.approx(-np.sqrt(np.pi), rel=1e-12)
    assert sol['delta_at'](np.array([np.pi / 2]))[0] == pytest.approx(-4.0, rel=1e-12)


# ---------------------------------------------------------------------------
# eq. (19) rectilinear SSB coefficient and Dugdale identities (zone solver)
# ---------------------------------------------------------------------------

def test_eq19_rectilinear_c1():
    for L in (0.05, 0.02, 0.01):
        assert cb.c1_rectilinear(L=L) == pytest.approx(np.sqrt(8 / np.pi), abs=2e-5)


def test_dugdale_delta_t():
    # delta_t = 8 sigma0 L/(pi E') for uniform closure, K_tip = 0 (eqs 11, 26)
    for L in (0.05, 0.02, 0.01):
        q, K0, Kt, info = ssb_zone_solve(L, s=0.0, sigma0=1.0, K0=0.0, N=48)
        assert -info['delta_t'] == pytest.approx(8 * L / np.pi, abs=2e-5)


# ---------------------------------------------------------------------------
# eq. (20): steady-state zone size with matrix toughness, two routes
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("g", [0.0, 0.1, 0.5, 1.0, 3.0])
def test_eq20_rectilinear_L0(g):
    # steady_state_rectilinear solves the zone integral equation numerically at
    # each trial L (no closed-form Dugdale coefficients), so agreement is at the
    # solver's accuracy (~1e-4..1e-3), not machine precision.
    L0_closed = cb.L0_rectilinear(g)
    L0_solver = cb.steady_state_rectilinear(g)
    assert L0_solver == pytest.approx(L0_closed, rel=5e-3)


def test_eq20_minus_form_and_planted_exponent():
    # (18)+(19) give the algebraically equivalent minus form
    # sqrt(L0) = [sqrt(1+g) - sqrt(g)] sqrt(pi/8); the printed plus-bracket
    # carries exponent -2.
    for g in (0.0, 0.1, 0.5, 1.0, 3.0):
        minus_form = (np.pi / 8.0) * (np.sqrt(1.0 + g) - np.sqrt(g))**2
        assert cb.L0_rectilinear(g) == pytest.approx(minus_form, rel=1e-12)
    # a planted +2 exponent must be rejected by the numerical route
    g = 0.1
    L0_wrong = (np.pi / 8.0) * (np.sqrt(1.0 + g) + np.sqrt(g))**2
    L0_numeric = cb.steady_state_rectilinear(g)
    assert abs(L0_numeric - L0_wrong) / L0_wrong > 0.2


def test_eq11_recovered_at_g0():
    assert cb.L0_rectilinear(0.0) == pytest.approx(np.pi / 8.0, rel=1e-12)


# ---------------------------------------------------------------------------
# Table 2 (hardening bridging) + eq. (26) identities
# ---------------------------------------------------------------------------

def test_eq26_identities():
    for beta in (0.09817, 0.39270, 0.78540, -0.07854, -0.11781, -0.3):
        c1, c2, c3, c4 = ssb_coefficients(beta, L=0.02, N=48)
        assert c2 == pytest.approx(np.sqrt(1 + 2 * beta * c1**2), rel=2e-3)
        assert c3 == pytest.approx((c2 - 1) / beta, rel=2e-3)
        assert c4 == pytest.approx(2 * c1, rel=2e-3)


def test_table2():
    # measured max errors of this solver vs paper: |dc1| = 0.0040, |dc2| = 0.0083
    rows = cb.reproduce_table2()
    for r in rows:
        assert r['c1'] == pytest.approx(r['c1_paper'], abs=0.005)
        assert r['c2_from_eq26'] == pytest.approx(r['c2_paper'], abs=0.009)


# ---------------------------------------------------------------------------
# Table 3 (softening bridging) and plateau eqs (27), (28)
# ---------------------------------------------------------------------------

def test_table3():
    # measured max errors of this solver vs paper: |dc1| = 0.0020, |d2G| = 0.0008
    rows = cb.reproduce_table3()
    for r in rows:
        assert r['c1'] == pytest.approx(r['c1_paper'], abs=0.003)
        assert r['two_G'] == pytest.approx(r['two_G_paper'], abs=0.002)


def test_softening_plateau_eqs27_28():
    lam, c1, gfrac = cb.softening_plateau()
    assert lam == pytest.approx(0.366, abs=2e-3)     # eq. (28)
    assert gfrac == pytest.approx(0.5, abs=5e-3)     # eq. (27)
    assert c1 == pytest.approx(1.169, abs=0.01)      # Table 3 last row


def test_table3_last_row_is_plateau():
    rows = cb.reproduce_table3()
    last = rows[-1]
    assert last['lam'] == pytest.approx(0.366, abs=3e-3)
    assert last['delta_t'] == pytest.approx(1.0, abs=5e-3)


# ---------------------------------------------------------------------------
# Fig. 8 R-curve kinks
# ---------------------------------------------------------------------------

def test_fig8_kinks():
    # rectilinear: straight R-curve of slope 8/pi up to the kink at pi/8, then 1
    x = np.linspace(0, 2, 200)
    y = cb.rcurve_rectilinear_x(x)
    x_rise = np.pi / 8.0 - 0.01
    assert cb.rcurve_rectilinear_x(np.array([x_rise]))[0] == \
        pytest.approx(x_rise / (np.pi / 8.0), rel=1e-12)
    assert cb.rcurve_rectilinear_x(np.array([np.pi / 8.0]))[0] == pytest.approx(1.0)
    # softening triangle (delta0 doubled to match Gamma): kink at 2*lam* = 0.732,
    # R-curve level 2G/(sigma0*delta0) = 1.000 at the plateau (Table 3 last row)
    lam, _, _ = cb.softening_plateau()
    assert 2 * lam == pytest.approx(0.732, abs=0.01)   # Fig. 8 softening kink
    xs, ys = cb.rcurve_softening_points()
    assert xs[-1] == pytest.approx(0.732, abs=0.01)
    assert ys[-1] == pytest.approx(1.0, abs=0.01)
    assert np.all(np.diff(ys) > 0)                       # still rising to the kink


# ---------------------------------------------------------------------------
# Strength: eqs (44), (47); spectral vs closed form (independent routes)
# ---------------------------------------------------------------------------

def test_strength_eq47_vs_spectral():
    # measured spectral error 0.06-0.65 % at N=700 (junction singularity)
    for alpha in (0.2, 0.5, 1.0, 1.4):
        s_closed = cb.strength_LSB(alpha)
        s_num = cb.strength_numerical(alpha, s=0.0, N=700, nL=30)
        assert s_num == pytest.approx(s_closed, rel=0.01)


def test_strength_SSB_LSB_asymptotics():
    # eq (47) approaches eq (44) from BELOW as alpha -> infinity and never
    # crosses it: gap ~9 % at alpha=0.7, ~3 % at alpha=2, ~0.65 % at alpha=10
    for alpha in (0.7, 2.0, 10.0):
        assert cb.strength_LSB(alpha) < cb.strength_SSB(alpha)
    assert cb.strength_LSB(0.7) == pytest.approx(0.6134, abs=0.002)
    assert cb.strength_SSB(0.7) == pytest.approx(0.6743, abs=0.002)
    gap = 1.0 - cb.strength_LSB(10.0) / cb.strength_SSB(10.0)
    assert gap == pytest.approx(0.0065, abs=0.001)


def test_strength_eq44_ssb():
    for alpha in (0.5, 1.0, 2.0):
        assert cb.strength_SSB(alpha) == pytest.approx(
            1 / np.sqrt(np.pi * alpha), rel=1e-12)


def test_strength_LSB_limit():
    # S/sigma0 -> 1 as alpha -> 0 (strength bounded by bridging strength)
    assert cb.strength_LSB(1e-4) == pytest.approx(1.0, abs=2e-4)


def test_strength_fig11_rectilinear_points():
    # descending (rectilinear) curve of Fig. 11, read from the figure;
    # closed form reads 0.910/0.699/0.528/0.455 -- within 0.011 of the reads
    fig11 = {0.2: (0.90, 0.03), 0.5: (0.70, 0.03), 1.0: (0.53, 0.03),
             1.4: (0.45, 0.03)}
    for a, (yref, tol) in fig11.items():
        assert cb.strength_LSB(a) == pytest.approx(yref, abs=tol)


def test_strength_fig11_softening_points():
    # ascending-curve = lower curve of Fig. 11 (softening triangle), read off;
    # spectral solver at test settings reads 0.642/0.446/0.361/0.313 --
    # within 0.017 of the figure reads (claimed figure accuracy +-0.03-0.05)
    fig11 = {0.2: (0.65, 0.05), 0.6: (0.44, 0.05), 1.0: (0.36, 0.05),
             1.4: (0.33, 0.05)}
    for a, (yref, tol) in fig11.items():
        s_num = cb.strength_numerical(a, s=-1.0, N=700, nL=30)
        assert s_num == pytest.approx(yref, abs=tol)


# ---------------------------------------------------------------------------
# Notch limits (Section 3.3)
# ---------------------------------------------------------------------------

def test_notch_limits():
    assert cb.notch_brittle_limit(0.0) == pytest.approx(1 / 3, rel=1e-12)
    aw = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    assert np.allclose(cb.notch_ductile_limit(aw), 1 - aw)


# ---------------------------------------------------------------------------
# Planted-error checks
# ---------------------------------------------------------------------------

def test_planted_error_wrong_c2_formula():
    """c2 = 1 + 2*beta*c1^2 (missing sqrt) must NOT match the solver c2."""
    beta = 0.39270
    c1, c2, c3, c4 = ssb_coefficients(beta, L=0.02, N=48)
    wrong = 1 + 2 * beta * c1**2
    assert not np.isclose(wrong, c2, rtol=0.02)


def test_planted_error_wrong_bridging_sign():
    """Flipping the spring sign changes c1 far outside the paper value."""
    # anti-hardening (s < 0 at beta = +0.098) must not reproduce Table 2 c1
    beta = 0.09817
    c1_wrong, _, _, _ = ssb_coefficients(-beta, L=0.02, N=48)
    paper_c1 = 1.670
    assert abs(c1_wrong - paper_c1) > 0.05


def test_planted_error_wrong_strength_formula():
    """Griffith formula (44) is not the LSB curve: must fail at small alpha."""
    assert not np.isclose(cb.strength_SSB(0.2), cb.strength_LSB(0.2), rtol=0.1)
