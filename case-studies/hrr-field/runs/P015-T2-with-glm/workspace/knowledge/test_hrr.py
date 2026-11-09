"""Tests for the Hutchinson (1968) HRR singular-field reproduction.

Run with:  python3 -m pytest test_hrr.py -v
"""
import inspect
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

import hrr

# ---------------------------------------------------------------------------
# Table 1 of the paper.  Values transcribed from the page image (p.21).
# Plane-stress I(n=13) is 2.87 in the image (paper.txt OCR wrongly says 2.37);
# plane-strain E_n = 1.23,1.37,1.48,1.52 (OCR of 1.43 and 1.37 is unreliable;
# the computed values 1.227,1.382,1.485,1.518 support these).
# ---------------------------------------------------------------------------
TABLE1 = {
    (3, False):  {"I": 3.86, "En": 2.16},
    (5, False):  {"I": 3.41, "En": 2.45},
    (9, False):  {"I": 3.03, "En": 2.67},
    (13, False): {"I": 2.87, "En": 2.76},
    (3, True):   {"I": 5.51, "En": 1.23},
    (5, True):   {"I": 5.01, "En": 1.37},
    (9, True):   {"I": 4.60, "En": 1.48},
    (13, True):  {"I": 4.40, "En": 1.52},
}

# Tolerances match the accuracy the unit claims: I reproduced to <=0.3%
# (worst observed ~0.27%) and E_n to <=0.85% (worst observed ~0.85%, plane
# strain n=5).  They are bounds on agreement with the printed table, not
# quadrature/ODE error estimates: the source of the historical table
# discrepancy is NOT established, since refining the present
# quadrature/ODE only changes E_n by ~1e-8, far below the 0.85% gap.
TOL_I = 0.005
TOL_EN = 0.010


@pytest.mark.parametrize("key", list(TABLE1))
def test_table1_I_matches_paper(key):
    n, pe = key
    sol = hrr.solve_angular(n, pe)
    I = hrr.integral_I(n, pe, sol)
    Ip = TABLE1[key]["I"]
    assert abs(I - Ip) / Ip < TOL_I, (n, pe, I, Ip)


@pytest.mark.parametrize("key", list(TABLE1))
def test_table1_En_matches_paper(key):
    n, pe = key
    sol = hrr.solve_angular(n, pe)
    I = hrr.integral_I(n, pe, sol)
    En = hrr.energy_coefficient_En(n, pe, sol, I)
    Ep = TABLE1[key]["En"]
    assert abs(En - Ep) / Ep < TOL_EN, (n, pe, En, Ep)


def test_energy_coefficient_refinement_value():
    """Pin the refinement numbers quoted in ENTRY (item 6 of the verification):
    plane-strain n=5 E_n is 1.3816415028 at 4001 points and 1.3816415137 at
    16001 points (change ~1.1e-8), far below the 0.85% table gap."""
    sol4 = hrr.solve_angular(5, True, n_theta=4001)
    E4 = hrr.energy_coefficient_En(5, True, sol4)
    sol16 = hrr.solve_angular(5, True, n_theta=16001)
    E16 = hrr.energy_coefficient_En(5, True, sol16)
    assert abs(E4 - 1.3816415028) < 1e-9, E4
    assert abs(E16 - 1.3816415137) < 1e-9, E16
    assert 0.0 < (E16 - E4) < 1e-7, (E4, E16)


def test_table1_amplitude_column():
    """(pi/I)^(1/(n+1)) (Table 1, second row of each block)."""
    paper = {
        (3, False): 0.949, (5, False): 0.987, (9, False): 1.004, (13, False): 1.006,
        (3, True): 0.869, (5, True): 0.925, (9, True): 0.963, (13, True): 0.976,
    }
    for (n, pe), val in paper.items():
        col = (np.pi / hrr.integral_I(n, pe)) ** (1.0 / (n + 1.0))
        assert abs(col - val) < 3e-3, (n, pe, col, val)


# ---------------------------------------------------------------------------
# Independent route 1: recover the exponent s = (2n+1)/(n+1) directly from the
# eigenvalue problem, without using the closed-form formula.
# ---------------------------------------------------------------------------
def _eigen_residual(sp, pe, n):
    s, p = sp
    if not (1.0 < s < 3.0):
        return [1e3, 1e3]
    try:
        sol = hrr._integrate(pe, n, p, s=s)
    except Exception:
        return [1e3, 1e3]
    return [sol.y[0, -1], sol.y[1, -1]]


@pytest.mark.parametrize("n", [3, 5, 9, 13])
@pytest.mark.parametrize("pe", [False, True])
def test_exponent_from_eigenvalue(n, pe):
    s_formula = hrr.singular_exponent(n)
    p0 = hrr.solve_angular(n, pe)["p"]
    r = least_squares(_eigen_residual, [s_formula * 1.01, p0], args=(pe, n),
                      xtol=1e-14, ftol=1e-14, gtol=1e-14)
    assert np.hypot(*r.fun) < 1e-9
    assert abs(r.x[0] - s_formula) < 1e-6, (n, pe, r.x[0], s_formula)


# ---------------------------------------------------------------------------
# Independent route 2: a known special solution -- the linear-elastic mode-I
# crack.  Its angular function is the biharmonic solution (eq. 40); the
# first-principles circular J integral must equal K_I^2/E.
#
# The integrand must include the derivatives of the polar basis vectors and
# the theta-derivatives of the displacement angular functions.  The earlier
# version of this test omitted utilde_r' and utilde_theta' and used the wrong
# cos(theta) coefficient (1.5 instead of a = s-1 = 0.5); those two mistakes
# cancelled and the test passed by accident (finding 2 of REVIEW.md).
# ---------------------------------------------------------------------------
def test_elastic_mode_I_J_equals_KI2_over_E_numeric():
    K_I, E, nu = np.sqrt(2.0 * np.pi), 1.0, 0.3
    th = np.linspace(-np.pi, np.pi, 200001)
    f, f1, f2, _ = hrr.elastic_angular_stress_function(th)
    sr = 1.5 * f + f2
    st = 0.75 * f
    srt = -0.5 * f1
    W = 0.5 * (sr ** 2 + st ** 2 - 2 * nu * sr * st + 2 * (1 + nu) * srt ** 2)
    ur = 0.5 * ((5 - 3 * nu) * np.cos(th / 2) - (1 + nu) * np.cos(1.5 * th))
    ut = 0.5 * ((nu - 7) * np.sin(th / 2) + (1 + nu) * np.sin(1.5 * th))
    urp = -0.25 * (5 - 3 * nu) * np.sin(th / 2) + 0.75 * (1 + nu) * np.sin(1.5 * th)
    utp = 0.25 * (nu - 7) * np.cos(th / 2) + 0.75 * (1 + nu) * np.cos(1.5 * th)
    a = 0.5  # s - 1 for s = 3/2
    integrand = (W * np.cos(th)
                 - a * np.cos(th) * (sr * ur + srt * ut)
                 - np.sin(th) * (sr * (ut - urp) - srt * (ur + utp)))
    J = (K_I ** 2 / (2.0 * np.pi)) * np.trapezoid(integrand, th) / E
    assert abs(J - K_I ** 2 / E) / (K_I ** 2 / E) < 1e-6, (J, K_I ** 2 / E)


def test_elastic_J_plane_stress_export_matches_numeric():
    """elastic_J_plane_stress (closed form) equals the first-principles
    calculus done by elastic_J_numeric, for several K_I and E."""
    for K_I, E in [(1.0, 1.0), (2.5, 3.0), (np.sqrt(2 * np.pi), 0.7)]:
        closed = hrr.elastic_J_plane_stress(K_I, E)
        numeric = hrr.elastic_J_numeric(K_I, E)
        assert abs(numeric - closed) / closed < 1e-6, (K_I, E, numeric, closed)


def test_elastic_angular_displacements_satisfy_constitutive_relations():
    """elastic_angular_displacements must follow from the angular stresses via
    the plane-stress constitutive law and the strain-displacement relations:
        ur = 2 (sr - nu st),
        ut = 2 (ur' - gamma),  gamma = 2 (1 + nu) srt   (engineering shear),
    checked at several theta against independently formed stresses."""
    nu = 0.3
    th = np.linspace(0.05, np.pi - 0.05, 13)
    f, fp, fpp, _ = hrr.elastic_angular_stress_function(th)
    sr = 1.5 * f + fpp
    st = 0.75 * f
    srt = -0.5 * fp
    ur, ut, urp, _ = hrr.elastic_angular_displacements(th, nu)
    assert np.max(np.abs(ur - 2.0 * (sr - nu * st))) < 1e-12
    gamma = 2.0 * (1.0 + nu) * srt
    assert np.max(np.abs(ut - 2.0 * (urp - gamma))) < 1e-12


def test_elastic_dimensionless_integral_is_2pi():
    """The explanatory equation of elastic_J_numeric carries the 1/(2*pi)
    normalisation: the dimensionless integral must be exactly 2*pi, so that
    J = (K_I^2/(2*pi*E)) * 2*pi = K_I^2/E (item 1 of the verification)."""
    th = np.linspace(-np.pi, np.pi, 200001)
    integrand, _ = hrr.elastic_J_integrand(th)
    N = np.trapezoid(integrand, th)
    assert abs(N - 2.0 * np.pi) < 1e-6, N
    # the export must be the integral of that same integrand (so the component
    # check below guards the pointwise integrand actually used by the export)
    K, E = 2.5, 3.0
    assert abs(hrr.elastic_J_numeric(K, E)
               - K ** 2 / (2 * np.pi * E) * N) < 1e-12
    # with the correct 1/(2*pi) prefactor, the unit field gives
    # J = K_I^2/E for K_I = sqrt(2*pi)
    assert abs(hrr.elastic_J_numeric(np.sqrt(2 * np.pi), 1.0) - 2 * np.pi) < 1e-5


def test_elastic_reduced_integrand_matches_cartesian_components():
    """Component check of the elastic export (V8 of the verification): the
    *actual returned* integrand of elastic_J_integrand must satisfy the
    pointwise identity
        integrand = W cos(theta) - [traction work from full Cartesian gradients]
    on the unit contour.  The traction work is evaluated independently from
    central differences of elastic_displacement_cartesian and the Cartesian
    stress tensor.  A wrong displacement-derivative expression in the returned
    integrand changes it pointwise even when the integrated J still cancels to
    K_I^2/E, so the earlier suite (which rebuilt the traction work from the
    intermediates and never read the returned value) let that regression pass.
    The intermediate-component assertion is retained below."""
    nu = 0.3
    th = np.linspace(-np.pi, np.pi, 9)
    integrand, T = hrr.elastic_J_integrand(th, nu)
    sr, st, srt = T["sr"], T["st"], T["srt"]
    ur, ut, urp, utp, a = T["ur"], T["ut"], T["urp"], T["utp"], T["a"]
    W = T["W"]
    cth, sth = np.cos(th), np.sin(th)
    reduced = (a * cth * (sr * ur + srt * ut)
               + sth * (sr * (ut - urp) - srt * (ur + utp)))
    x, y = cth, sth            # r0 = 1
    hh = 1e-6

    def dux_dx(px, py):
        return (hrr.elastic_displacement_cartesian(px + hh, py, nu, 1.0)[0]
                - hrr.elastic_displacement_cartesian(px - hh, py, nu, 1.0)[0]) / (2 * hh)

    def duy_dx(px, py):
        return (hrr.elastic_displacement_cartesian(px + hh, py, nu, 1.0)[1]
                - hrr.elastic_displacement_cartesian(px - hh, py, nu, 1.0)[1]) / (2 * hh)

    sxx = sr * cth ** 2 + st * sth ** 2 - 2 * srt * sth * cth
    sxy = (sr - st) * sth * cth + srt * (cth ** 2 - sth ** 2)
    syy = sr * sth ** 2 + st * cth ** 2 + 2 * srt * sth * cth
    tx = sxx * cth + sxy * sth
    ty = sxy * cth + syy * sth
    cart = np.array([tx[i] * dux_dx(x[i], y[i]) + ty[i] * duy_dx(x[i], y[i])
                     for i in range(len(th))])
    # intermediate components still agree with the Cartesian traction work
    assert np.max(np.abs(reduced - cart)) < 1e-8, np.max(np.abs(reduced - cart))
    # and the *returned* integrand is W cos(theta) minus that same traction work
    resid = np.max(np.abs(integrand - (W * cth - cart)))
    assert resid < 1e-8, resid


def test_elastic_displacement_cartesian_modulus_scaling():
    """elastic_displacement_cartesian must carry the 1/E modulus factor
    (item 5 of the verification): for the point (x,y)=(1,0) and nu=0.3 the
    E=1 displacement is u=(1.4, 0), so E=3 must give (1.4/3, 0).  A helper that
    omits the division by E passes every test that only calls it with E=1."""
    x, y, nu = 1.0, 0.0, 0.3
    ux1, uy1 = hrr.elastic_displacement_cartesian(x, y, nu, 1.0)
    ux3, uy3 = hrr.elastic_displacement_cartesian(x, y, nu, 3.0)
    assert abs(ux1 - 1.4) < 1e-12, ux1
    assert abs(uy1) < 1e-12, uy1
    assert abs(ux3 - 1.4 / 3.0) < 1e-12, ux3
    assert abs(uy3) < 1e-12, uy3
    # a second, off-axis point, checked against an independent scaling
    x2, y2 = 0.6, 0.8
    u1 = np.array(hrr.elastic_displacement_cartesian(x2, y2, nu, 1.0))
    uE = np.array(hrr.elastic_displacement_cartesian(x2, y2, nu, 2.5))
    assert np.max(np.abs(uE - u1 / 2.5)) < 1e-14, (u1, uE)


def test_elastic_angular_stress_function_known_values():
    """eq. (40): f(0)=4/3, f(pi)=f'(pi)=f''(pi)=0, f''(0)=-1, f'''(pi)=-1."""
    f, fp, fpp, fppp = hrr.elastic_angular_stress_function(
        np.array([0.0, np.pi]))
    assert abs(f[0] - 4.0 / 3.0) < 1e-12
    assert abs(fp[0]) < 1e-12
    assert abs(fpp[0] + 1.0) < 1e-12
    assert abs(f[1]) < 1e-12
    assert abs(fp[1]) < 1e-12
    assert abs(fppp[1] + 1.0) < 1e-12
    # traction-free crack faces: eq. (14) with s = 3/2 gives
    # sr = s f + f'' and srt = (1-s) f' both zero at theta = pi.
    s = 1.5
    assert abs(s * f[1] + fpp[1]) < 1e-12
    assert abs((1 - s) * fp[1]) < 1e-12


# ---------------------------------------------------------------------------
# Route 3: the tabulated energy identity
#   I * E_n = (n/(n+1)) * pi * int_{-pi}^{pi} se_tilde^(n+1) dtheta
# This does NOT independently validate I: E_n is defined from the same I, so I
# cancels algebraically.  What it validates is (a) the internal consistency of
# the E_n definition and (b) the angular energy shape through its product with
# the paper's I*E_n (the shape part survives a rescaling of I).  It is not a
# second independent J route (finding 12 of REVIEW.md).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("key", list(TABLE1))
def test_energy_shape_identity(key):
    n, pe = key
    sol = hrr.solve_angular(n, pe)
    shape = 2.0 * np.trapezoid(sol["se"] ** (n + 1.0), sol["theta"])
    lhs = hrr.integral_I(n, pe, sol) * hrr.energy_coefficient_En(n, pe, sol)
    rhs = (n / (n + 1.0)) * np.pi * shape
    assert abs(lhs - rhs) / rhs < 1e-10
    got = lhs
    paper = TABLE1[key]["I"] * TABLE1[key]["En"]
    assert abs(got - paper) / paper < 0.02


# ---------------------------------------------------------------------------
# Independent route for I: evaluate Rice's J integral from the similarity field
# with *Cartesian* finite differences (so the derivatives of the polar basis are
# included automatically), and compare with the angular quadrature of eq. (24).
# This is the route that refutes the former "I = 5.65" discrepancy (findings
# 10/N2): the earlier direct route had dropped the angular derivatives.
# ---------------------------------------------------------------------------
def _cartesian_J_from_solution(n, pe, sol, r0=1.0, hh=1e-5):
    s, th = sol["s"], sol["theta"]
    f, fp, fpp, fppp = sol["f"], sol["fp"], sol["fpp"], sol["fppp"]
    se, sr, st, srt = sol["se"], sol["sr"], sol["st"], sol["srt"]
    M = hrr.effective_stress_coeff(n, pe)
    z = np.stack([f, fp, fpp], -1)
    dQ = 2.0 * (z @ M.T)
    dQv = np.einsum('ij,ij->i', dQ, np.stack([fp, fpp, fppp], -1))
    rho = se ** (n - 1.0)
    rhop = (n - 1.0) * rho * (dQv / (2.0 * se ** 2))
    m = n * (s - 2.0)
    b = m + 1.0
    if pe:
        c1 = 0.75 * (sr - st)
        c1p = 0.75 * (s * (2.0 - s) * fp + fppp)
    else:
        c1 = sr - 0.5 * st
        c1p = (s * (3.0 - s) / 2.0) * fp + fppp
    ur = rho * c1 / b
    urp = (rhop * c1 + rho * c1p) / b
    ut = (3.0 * rho * srt - urp) / m
    if pe:
        utp = rho * (0.75 * (st - sr)) - ur
    else:
        utp = rho * (s * (s - 1.5) * f - fpp / 2.0) - ur
    a = n * (s - 2.0) + 1.0

    def disp(x, y):
        r = np.hypot(x, y)
        t = np.arctan2(y, x)
        ta = np.abs(t)
        urv = r ** a * np.interp(ta, th, ur)
        utv = r ** a * np.interp(ta, th, ut) * np.sign(t)
        ct, stt = x / r, y / r
        return urv * ct - utv * stt, urv * stt + utv * ct

    def stress(x, y):
        r = np.hypot(x, y)
        t = np.arctan2(y, x)
        ta = np.abs(t)
        fac = r ** (s - 2.0)
        return (fac * np.interp(ta, th, sr),
                fac * np.interp(ta, th, srt) * np.sign(t))

    ct, stt = np.cos(th), np.sin(th)
    x, y = r0 * ct, r0 * stt
    uxp, uyp = disp(x + hh, y)
    uxm, uym = disp(x - hh, y)
    dux_dx = (uxp - uxm) / (2.0 * hh)
    duy_dx = (uyp - uym) / (2.0 * hh)
    sr_, srt_ = stress(x, y)
    tx = sr_ * ct - srt_ * stt
    ty = sr_ * stt + srt_ * ct
    W = (n / (n + 1.0)) * r0 ** ((n + 1.0) * (s - 2.0)) * se ** (n + 1.0)
    integrand = (W * ct - (tx * dux_dx + ty * duy_dx)) * r0
    return 2.0 * np.trapezoid(integrand, th)


@pytest.mark.parametrize("n,pe", [(3, False), (3, True), (9, True)])
@pytest.mark.parametrize("r0", [0.2, 1.0, 3.0])
def test_integral_I_matches_cartesian_J(n, pe, r0):
    sol = hrr.solve_angular(n, pe)
    I = hrr.integral_I(n, pe, sol)
    J = _cartesian_J_from_solution(n, pe, sol, r0=r0)
    assert abs(J - I) / I < 1e-5, (n, pe, r0, I, J)


# ---------------------------------------------------------------------------
# Exported functions required to be callable (finding 11 of REVIEW.md).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("pe", [False, True])
def test_effective_stress_coeff_matches_direct_formula(pe):
    """effective_stress_coeff must reproduce se^2 from the stress components:
    plane stress : se^2 = sr^2 + st^2 - sr st + 3 srt^2
    plane strain : se^2 = (3/4)(sr-st)^2 + 3 srt^2.
    """
    n = 3
    sol = hrr.solve_angular(n, pe)
    M = hrr.effective_stress_coeff(n, pe)
    z = np.stack([sol["f"], sol["fp"], sol["fpp"]], -1)
    se_M2 = np.einsum('ij,jk,ik->i', z, M, z)
    if pe:
        direct = 0.75 * (sol["sr"] - sol["st"]) ** 2 + 3.0 * sol["srt"] ** 2
    else:
        direct = (sol["sr"] ** 2 + sol["st"] ** 2 - sol["sr"] * sol["st"]
                  + 3.0 * sol["srt"] ** 2)
    assert np.max(np.abs(se_M2 - direct)) < 1e-12
    # and it must be the matrix that normalises max se_tilde to 1
    assert abs(np.sqrt(np.max(se_M2)) - 1.0) < 1e-12


@pytest.mark.parametrize("pe", [False, True])
def test_amplitude_K_satisfies_J_balance(pe):
    """amplitude_K must satisfy the small-scale-yielding J balance

        alpha K^(n+1) I = pi sigma_inf^2            (plane stress)
        alpha K^(n+1) I = pi (1-nu^2) sigma_inf^2   (plane strain)

    using the *paper's* I from Table 1 (i.e. an independent value, not the one
    the function used).  A corrupted amplitude_K cannot satisfy this.
    """
    n, alpha, sigma_inf, nu = 3, 2.0, 1.5, 0.3
    I_paper = TABLE1[(n, pe)]["I"]
    K = hrr.amplitude_K(n, pe, alpha=alpha, sigma_inf=sigma_inf, nu=nu,
                        I=I_paper)
    J_local = alpha * K ** (n + 1.0) * I_paper
    J_far = np.pi * sigma_inf ** 2 * (1.0 if not pe else (1.0 - nu ** 2))
    assert abs(J_local - J_far) / J_far < 1e-12, (pe, K, J_local, J_far)
    # explicit paper formula cross-check
    amp = 1.0 if not pe else (1.0 - nu ** 2)
    K_ref = ((amp / alpha) ** (1.0 / (n + 1.0))
             * (np.pi / I_paper) ** (1.0 / (n + 1.0))
             * sigma_inf ** (2.0 / (n + 1.0)))
    assert abs(K - K_ref) / K_ref < 1e-14


@pytest.mark.parametrize("pe", [False, True])
def test_amplitude_K_documented_reconstruction(pe):
    """Call amplitude_K exactly as corrected ENTRY.md documents it and check the
    documented reconstruction K = (J_dimless/(alpha*I))^(1/(n+1)) for BOTH
    constraints (V4 of the verification).  J_dimless is the *outer* elastic J,
    pi sigma_inf^2 (plane stress) or pi (1-nu^2) sigma_inf^2 (plane strain); the
    1-nu^2 factor is already inside J and must not be divided out again when
    inverting for K.  The round-2 documentation's extra 1/(1-nu^2) gives a K
    that fails this J balance."""
    n, alpha, sigma_inf, nu = 3, 2.0, 0.2, 0.3
    I_paper = TABLE1[(n, pe)]["I"]
    K = hrr.amplitude_K(n, plane_strain=pe, alpha=alpha,
                        sigma_inf=sigma_inf, nu=nu, I=I_paper)
    J_outer = np.pi * sigma_inf ** 2 * (1.0 if not pe else (1.0 - nu ** 2))
    K_recon = (J_outer / (alpha * I_paper)) ** (1.0 / (n + 1.0))
    assert abs(K - K_recon) < 1e-14, (pe, K, K_recon)
    # the reconstruction really is a J balance: alpha K^(n+1) I = J_outer
    assert abs(alpha * K ** (n + 1.0) * I_paper - J_outer) < 1e-14
    if pe:
        K_bad = (J_outer / ((1.0 - nu ** 2) * alpha * I_paper)) ** (1.0 / (n + 1.0))
        assert abs(K - K_bad) / K > 0.02, (K, K_bad)
        # the review's worked numbers with the paper's I = 5.51:
        # J = pi(1-nu^2)sigma_inf^2 = 0.11435397..., K = 0.31916658...
        assert abs(J_outer - 0.1143539726) < 1e-9, J_outer
        assert abs(K - 0.31916658) < 1e-7, K


def test_documented_signatures_and_usage_statements():
    """Every documented call form and usage statement in ENTRY.md's
    "Functions to call" / "Optional arguments" text must match the real API
    (G1 of the round-3 verification: ENTRY claimed `integral_I` accepts `I`,
    which it does not).

    The test (a) calls each of the 12 documented signatures exactly as written,
    (b) checks the argument-binding claims against reality, and (c) reads
    ENTRY.md, so it fails if the false instruction is reintroduced.
    """
    entry = Path(__file__).with_name("ENTRY.md").read_text()
    # the corrected statement is present and the round-3 false one is gone
    assert ("`integral_I`, `energy_coefficient_En` and `amplitude_K` all accept `I`"
            not in entry)
    assert ("Only `energy_coefficient_En` and `amplitude_K` accept a supplied `I`"
            in entry)

    # documented signature block vs the real API
    sigs = {name: inspect.signature(getattr(hrr, name)) for name in hrr.__all__}
    assert list(sigs["integral_I"].parameters) == ["n", "plane_strain", "sol"]
    assert list(sigs["energy_coefficient_En"].parameters) == [
        "n", "plane_strain", "sol", "I"]
    assert list(sigs["amplitude_K"].parameters) == [
        "n", "plane_strain", "alpha", "sigma_inf", "nu", "I"]

    # every documented signature called exactly as written
    n, pe, nu = 3, False, 0.3
    th = np.array([0.0, 0.5, 1.0])
    s = hrr.singular_exponent(n)
    M = hrr.effective_stress_coeff(n, plane_strain=pe)
    sol = hrr.solve_angular(n, plane_strain=pe, s=None, normalize=True)
    I = hrr.integral_I(n, plane_strain=pe, sol=sol)
    En = hrr.energy_coefficient_En(n, plane_strain=pe, sol=sol, I=I)
    K = hrr.amplitude_K(n, plane_strain=pe, alpha=2.0, sigma_inf=0.2, nu=nu, I=I)
    f = hrr.elastic_angular_stress_function(th)
    u = hrr.elastic_angular_displacements(th, nu=nu)
    integrand, terms = hrr.elastic_J_integrand(th, nu=nu)
    ux = hrr.elastic_displacement_cartesian(1.0, 0.5, nu=nu, E=1.0)
    J_ps = hrr.elastic_J_plane_stress(1.0, 1.0)
    J_num = hrr.elastic_J_numeric(1.0, E=1.0, nu=nu)

    # the documented I-accepting calls bind and give finite positive results
    assert M.shape == (3, 3) and s > 1.0
    assert I > 0.0 and En > 0.0 and K > 0.0
    assert ux[0] > 0.0 and J_ps > 0.0 and J_num > 0.0
    assert np.isfinite(f).all() and np.isfinite(u).all()
    assert np.isfinite(integrand).all() and terms["W"].shape == th.shape

    # `integral_I` accepts `sol` (above) but, as documented, does not accept `I`
    with pytest.raises(TypeError):
        hrr.integral_I(n, plane_strain=pe, I=I)
    # `amplitude_K`, as documented, does not accept the `sol` keyword
    with pytest.raises(TypeError):
        hrr.amplitude_K(n, plane_strain=pe, sol=sol)


@pytest.mark.parametrize("pe", [False, True])
def test_integral_I_scales_as_field_to_n_plus_1(pe):
    """Normalisation scope (V4): with the raw field = normalised / c_norm, the
    dimensionless I scales as the angular field scale to the power n+1, i.e.
    I_raw = I_norm / c_norm^(n+1).  The standard table values hold only for the
    normalised field (max se = 1).  For n=3 plane stress this is
    3.854966353 (normalised) vs 7.643514258 (raw), c_norm = 1/1.1866376456."""
    n = 3
    nrm = hrr.solve_angular(n, pe, normalize=True)
    raw = hrr.solve_angular(n, pe, normalize=False)
    I_norm = hrr.integral_I(n, pe, nrm)
    I_raw = hrr.integral_I(n, pe, raw)
    c = nrm["c"]
    assert abs(I_raw - I_norm / c ** (n + 1.0)) < 1e-9, (pe, I_norm, I_raw, c)
    if not pe:
        # concrete values quoted by the verification; here c_norm < 1
        assert abs(I_norm - 3.854966353) < 1e-6, I_norm
        assert abs(I_raw - 7.643514258) < 1e-6, I_raw
        assert I_raw > I_norm        # raw field is larger, c_norm < 1
    # the raw and normalised I differ (the normalisation matters), and only the
    # normalised field gives the standard table value
    assert abs(I_raw - I_norm) / I_norm > 1e-3, (pe, I_norm, I_raw)


@pytest.mark.parametrize("pe", [False, True])
def test_amplitude_K_default_I_satisfies_J_balance(pe):
    """The *default* path (no I argument) of amplitude_K must also satisfy the
    small-scale-yielding J balance, with I computed internally.  This closes the
    default-I coverage gap (item 7 of the verification): an amplitude_K that
    swaps the constraint only for its internally computed I survives every test
    that always passes I explicitly."""
    n, alpha, sigma_inf, nu = 3, 2.0, 0.2, 0.3
    K = hrr.amplitude_K(n, pe, alpha=alpha, sigma_inf=sigma_inf, nu=nu)
    I_default = hrr.integral_I(n, pe)
    J_local = alpha * K ** (n + 1.0) * I_default
    J_far = np.pi * sigma_inf ** 2 * (1.0 if not pe else (1.0 - nu ** 2))
    assert abs(J_local - J_far) / J_far < 1e-12, (pe, K, I_default, J_local, J_far)


@pytest.mark.parametrize("pe", [False, True])
def test_solve_angular_normalize_false_returns_raw_field(pe):
    """normalize=False must leave c = 1 and return the un-normalised field:
    raw = normalised / c_norm and max(se_raw) = 1 / c_norm.  An implementation
    that ignores normalize=False survives tests that only look at the
    normalised field (item 7 of the verification)."""
    n = 3
    raw = hrr.solve_angular(n, pe, normalize=False)
    nrm = hrr.solve_angular(n, pe, normalize=True)
    assert raw["c"] == 1.0
    assert abs(raw["se"].max() - 1.0 / nrm["c"]) < 1e-12
    if not pe:
        # the concrete raw value quoted by the verification for n=3, plane stress
        assert abs(raw["se"].max() - 1.1866376456455583) < 1e-12
    # the shooting value p is unchanged by the normalisation, and the raw
    # angular functions are the normalised ones divided by c (the c/p scaling)
    assert abs(raw["p"] - nrm["p"]) < 1e-12
    assert np.max(np.abs(raw["f"] - nrm["f"] / nrm["c"])) < 1e-12
    assert np.max(np.abs(raw["sr"] - nrm["sr"] / nrm["c"])) < 1e-12
    # the normalised field really has max se_tilde = 1
    assert abs(nrm["se"].max() - 1.0) < 1e-12


@pytest.mark.parametrize("n", [3, 5, 9, 13])
def test_j_scaling_exponent_at_eigenvalue(n):
    """At the eigenvalue (n+1)(s-2) = -1, so W ~ 1/r (not r-independent); the
    extra factor r from the contour makes the J integrand r-independent
    (finding 4 of REVIEW.md)."""
    s = hrr.singular_exponent(n)
    assert abs((n + 1.0) * (s - 2.0) + 1.0) < 1e-14  # J integrand r^0
    assert (n + 1.0) * (s - 2.0) < 0.0               # W ~ r^(negative)
    assert abs((n + 1.0) * (s - 2.0)) > 0.5          # exponent is -1, not 0


# ---------------------------------------------------------------------------
# Planted-error checks: a deliberately wrong exponent / wrong effective stress
# must break the reproduction.
# ---------------------------------------------------------------------------
def test_planted_wrong_exponent_fails():
    """Using s = 2n/(n+1) instead of (2n+1)/(n+1) must be detected: either no
    physical root exists, or the resulting I is far outside Table 1."""
    orig = hrr.singular_exponent
    detected = False
    try:
        hrr.singular_exponent = lambda n: 2.0 * n / (n + 1.0)  # wrong
        try:
            I_wrong = hrr.integral_I(3, False)
            detected = abs(I_wrong - 3.86) / 3.86 > TOL_I
        except RuntimeError:
            detected = True
    finally:
        hrr.singular_exponent = orig
    assert detected


def test_planted_wrong_effective_stress_fails():
    """Using the plane-stress effective stress in the plane-strain ODE is wrong."""
    orig = hrr._quadratic_coeff

    def wrong_coeff(plane_strain, n, s):
        return orig(False, n, s)  # always plane stress

    try:
        hrr._quadratic_coeff = wrong_coeff
        I_wrong = hrr.integral_I(3, True)
    finally:
        hrr._quadratic_coeff = orig
    assert abs(I_wrong - 5.51) / 5.51 > TOL_I, I_wrong


def test_correct_value_passes_same_tolerance():
    """Sanity: the un-planted computation passes the same tolerance test."""
    assert abs(hrr.integral_I(3, False) - 3.86) / 3.86 < TOL_I
    assert abs(hrr.integral_I(3, True) - 5.51) / 5.51 < TOL_I
