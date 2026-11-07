#!/usr/bin/env python3
"""
Reference reproduction for PROBLEM.md ("Energy release for periodic cracks in an elastic film"),
limited to the two array geometries (Case A simultaneous array; Case B sequential midway cracking).

Single deterministic command:  python3 reproduce.py

Routes:
  (1) Analytic closed forms derived from the problem's stated incremental energy and
      total-traction boundary conditions (checked symbolically with SymPy if available).
  (2) Numerical strip boundary-value solution of v'' = v, v'(0) = -1 (total traction-free crack
      face: Delta_sigma22 = -sigma0), v(a) = 0 (midline antisymmetry), by collocation
      (scipy.integrate.solve_bvp).
  (3) Independent second-order finite-difference Newton solve of the same strip BVP.
The deciding energies are obtained by direct numerical quadrature of the stated dimensionless
functional I(S) = Int_0^{2a} [v' + (v'**2 + v**2)/2] deta = -tanh(a)  (a = S/(2*ell)),
using ONLY numerical v (never the closed form). Then
  gamma_A(S) = -I(S),   gamma_B = 2*gamma_A_num(H/ell) - gamma_A_num(2*H/ell).

Controls (deliberately WRONG, preserved with their actual failure outputs; NOT production answers):
  C1: omit the foundation energy (k/2) u^2 from the stated functional.
  C2: impose the INCREMENTAL traction Delta_sigma22 = 0 on the crack faces instead of the total
      traction sigma22 = 0 (i.e. v'(0) = 0 instead of -1).
  C3: confuse the initial and final strip spacings (swap the ahead/behind states of Case B).

Outputs: full stdout report, reference.json in this folder. Exit code 0 iff all production checks
pass AND all controls are confirmed to fail by a clear margin.
"""
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_bvp

HERE = os.path.dirname(os.path.abspath(__file__))
NU = 0.3
ONE_MINUS_NU2 = 1.0 - NU * NU
SAMPLES_ELL_OVER_H = [0.5, 1.0, 1.5, 2.0]

# ----------------------------------------------------------------------------
# Route 1: analytic closed forms (derived in REFERENCE.md Sections 1-5)
# ----------------------------------------------------------------------------

def gamma_A_analytic(S_over_ell):
    """G_A(S)*Ebar/(sigma0^2 ell) for an array of semi-infinite cracks with spacing S."""
    return math.tanh(S_over_ell / 2.0)

def gamma_B_analytic(H_over_ell):
    """G_B(H)*Ebar/(sigma0^2 ell) for sequential midway cracking between spacing-2H cracks."""
    return 2.0 * math.tanh(H_over_ell / 2.0) - math.tanh(H_over_ell)

def kI_from_gamma(g):
    """KI/(sigma0 sqrt(ell)) from gamma via pure mode I and G = KI^2/E, Ebar = E/(1-nu^2)."""
    return math.sqrt(ONE_MINUS_NU2 * g)

# ----------------------------------------------------------------------------
# Route 2: collocation strip BVP
# ----------------------------------------------------------------------------

def strip_bvp_collocation(a, tol=1e-11, n0=64):
    """Solve v'' = v on [0,a], v'(0) = -1, v(a) = 0. Returns solve_bvp solution object."""

    def fun(x, y):
        return np.vstack((y[1], y[0]))

    def bc(ya, yb):
        return np.array([ya[1] + 1.0, yb[0]])

    x = np.linspace(0.0, a, n0)
    th = math.tanh(a)
    y0 = np.vstack((th * (1.0 - x / a) ** 2, -2.0 * th * (1.0 - x / a) / a))  # generic start
    sol = solve_bvp(fun, bc, x, y0, tol=tol, max_nodes=1000000)
    if not sol.success:
        raise RuntimeError("collocation BVP failed: " + sol.message)
    return sol

def v_exact(eta, a):
    return np.sinh(a - eta) / np.cosh(a)

# ----------------------------------------------------------------------------
# Route 3: independent finite-difference strip BVP (2nd order, tridiagonal)
# ----------------------------------------------------------------------------

def strip_bvp_fd(a, N):
    """Return (eta, v) for the same BVP by a ghost-node finite-difference (linear) solve.
    Nodes i=0..N on [0,a]; v'(0) = -1 via ghost v_{-1} = v_1 + 2d; v(a) = 0 as the last row."""
    d = a / N
    eta = np.linspace(0.0, a, N + 1)
    diag = np.zeros(N + 1)
    sup = np.zeros(N)   # sup[i] = A[i, i+1], i = 0..N-1
    sub = np.zeros(N)   # sub[i] = A[i+1, i], i = 0..N-1
    rhs = np.zeros(N + 1)
    diag[0] = -(2.0 + d * d)      # (2 v1 - (2 + d^2) v0)/d^2 = 0 with ghost eliminated, rhs -2 d
    sup[0] = 2.0
    rhs[0] = -2.0 * d
    for i in range(1, N):
        diag[i] = -(2.0 + d * d)  # v_{i+1} - (2 + d^2) v_i + v_{i-1} = 0
        sub[i - 1] = 1.0
        sup[i] = 1.0
    diag[N] = 1.0                 # v_N = 0
    sub[N - 1] = 0.0
    from scipy.linalg import solve_banded

    ab = np.zeros((3, N + 1))
    ab[0, 1:] = sup
    ab[1, :] = diag
    ab[2, :-1] = sub
    v = solve_banded((1, 1), ab, rhs)
    v[N] = 0.0
    return eta, v

# ----------------------------------------------------------------------------
# Energy quadrature of the stated functional using ONLY a numerical v
# ----------------------------------------------------------------------------

def integrate_energy(eta_half, v_half, vp_half):
    """I = Int_0^{2a}[v' + (v'^2 + v^2)/2] deta using the half-strip solution and
    antisymmetry v(2a-x) = -v(x):  I = -2 v(0) + Int_0^a (v'^2 + v^2) deta."""
    dens = vp_half**2 + v_half**2
    # composite Simpson
    integ = np.sum(simpson_weights(eta_half) * dens)
    return -2.0 * v_half[0] + integ

def simpson_weights(x):
    n = len(x) - 1
    if n % 2 == 1:
        raise ValueError("Simpson needs an even number of panels")
    w = np.ones(n + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    h = x[1] - x[0]
    return w * h / 3.0

def gamma_num_from_solution(sol, S_over_ell, npts=4001):
    a = S_over_ell / 2.0
    eta = np.linspace(0.0, a, npts)
    v = sol.sol(eta)[0]
    vp = sol.sol(eta, 1)[0]
    I = integrate_energy(eta, v, vp)
    return -I, I

# ----------------------------------------------------------------------------
# Symbolic verification (Route 0): exact strip solution satisfies ODE/BCs and the
# stated functional integrates to -h sigma0^2 ell tanh(a)/Ebar per period.
# ----------------------------------------------------------------------------

def sympy_check():
    """Symbolically verify (i) the exact strip field solves the dimensional model with k = Ebar h/ell^2,
    (ii) total-traction face and midline conditions, (iii) the dimensional reduction of the stated
    functional to the dimensionless integrand, and (iv) the period integral -h sigma0^2 ell tanh(a)/Ebar."""
    try:
        import sympy as sp
    except Exception as e:  # pragma: no cover
        return ("SKIPPED (sympy unavailable: %s)" % e, None)
    eta, a, s0, Eb, h, ell = sp.symbols("eta a sigma0 Ebar h ell", positive=True)
    v = sp.sinh(a - eta) / sp.cosh(a)                     # dimensionless strip field on [0, 2a]
    # (i) dimensionless BVP and conditions
    ode_ok = sp.simplify(sp.diff(v, eta, 2) - v) == 0
    face_ok = sp.simplify(sp.diff(v, eta).subs(eta, 0) + 1) == 0   # Ebar u'(0) = -sigma0 (TOTAL traction)
    mid_ok = sp.simplify(v.subs(eta, a)) == 0
    # (ii) dimensional reduction: u = (sigma0*ell/Ebar) v, x2 = ell*eta, k = Ebar*h/ell**2
    #   equilibrium Ebar u'' = k u  and  integrand*dx2 = ell*(h sigma0^2/Ebar)*[v' + (v'^2+v^2)/2]*deta
    k = Eb * h / ell**2
    u = s0 * ell / Eb * v                                    # dimensional wake field u(x2 = ell*eta)
    eq_ok = sp.simplify(Eb * sp.diff(u, eta, 2) / ell**2 - k * u / h) == 0  # h*Ebar u_{,x2x2} = k u
    lhs = (h * s0 * sp.diff(u, eta) / ell
           + h / 2 * Eb * (sp.diff(u, eta) / ell) ** 2 + k / 2 * u**2) * ell   # integrand * dx2
    rhs = ell * h * s0**2 / Eb * (sp.diff(v, eta) + (sp.diff(v, eta) ** 2 + v**2) / 2)
    red_ok = sp.simplify(lhs - rhs) == 0
    # (iii) period integral of the stated functional: I = Int_0^{2a}[v' + (v'^2+v^2)/2] = -tanh(a)
    I = sp.integrate(sp.diff(v, eta) + (sp.diff(v, eta) ** 2 + v**2) / 2, (eta, 0, 2 * a))
    I_ok = sp.simplify(I + sp.tanh(a)) == 0
    ok = bool(ode_ok and face_ok and mid_ok and eq_ok and red_ok and I_ok)
    return (
        "symbolic checks: ODE=%s face=%s midline=%s equilibrium(k=Ebar*h/ell^2)=%s "
        "functional_reduction=%s period_integral=-tanh(a):%s" % (ode_ok, face_ok, mid_ok, eq_ok, red_ok, I_ok),
        ok,
    )

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("=" * 100)
    say("Reference reproduction: periodic cracks in a prestressed film on a Winkler foundation")
    say("model per PROBLEM.md; source: Xia & Hutchinson, JMPS 48 (2000) 1107-1131, Sec. 4.2, Eqs. 28, 30")
    say("nu = %.4f  (1-nu^2 = %.10f);  ell = sqrt(Ebar h / k);  Ebar = E/(1-nu^2)" % (NU, ONE_MINUS_NU2))
    say("=" * 100)

    # ---- symbolic check ------------------------------------------------------
    msg, ok = sympy_check()
    say("[sym] %s" % msg)
    say("")

    # ---- BVP validation and convergence -------------------------------------
    say("[BVP] strip problem v'' = v, v'(0) = -1 (TOTAL traction-free face), v(a) = 0 (midline)")
    say("[BVP] collocation (solve_bvp) max|v_num - sinh(a-eta)/cosh(a)| and Simpson-integrated")
    say("[BVP] energy residual I_num + tanh(a), plus finite-difference route error:")
    test_S = [1.0, 2.0, 2.0 / 3.0, 0.5, 4.0, 3.0]
    max_bvp_err = 0.0
    max_I_err = 0.0
    conv_rows = []
    for S_over_ell in sorted(set(test_S)):
        a = S_over_ell / 2.0
        sol = strip_bvp_collocation(a)
        eta = np.linspace(0.0, a, 20001)
        err = float(np.max(np.abs(sol.sol(eta)[0] - v_exact(eta, a))))
        I_num, I = gamma_num_from_solution(sol, S_over_ell)
        I_err = abs(I_num - math.tanh(a))
        max_bvp_err = max(max_bvp_err, err)
        max_I_err = max(max_I_err, I_err)
        # finite-difference route at several N (convergence of gamma)
        errs_fd = []
        for N in (80, 160, 320, 640, 1280):
            etaf, vf = strip_bvp_fd(a, N)
            vpf = np.gradient(vf, etaf, edge_order=2)
            I_fd = integrate_energy(etaf, vf, vpf)
            errs_fd.append(abs(-I_fd - math.tanh(a)))
        conv_rows.append((S_over_ell, err, I_err, errs_fd))
        say("  S/ell=%8.4f  colloc_max_err=%.3e  gamma_num_err=%.3e  fd_gamma_err(N=80..1280)=%s"
            % (S_over_ell, err, I_err, " ".join("%.2e" % e for e in errs_fd)))
    say("")
    tol_bvp = 1e-8
    assert max_bvp_err < tol_bvp, "collocation BVP error too large: %g" % max_bvp_err
    assert max_I_err < 1e-9, "energy quadrature residual too large: %g" % max_I_err
    fd_last = max(r[3][-1] for r in conv_rows)
    say("[BVP] checks passed: collocation max error %.3e < %g; energy residual %.3e; FD(1280) %.3e"
        % (max_bvp_err, tol_bvp, max_I_err, fd_last))
    say("")

    # ---- production table ----------------------------------------------------
    say("=" * 100)
    say("[MAIN] normalized results at nu = 0.3 (ell/H sampling locations, not material data)")
    say("%9s %9s %20s %20s %20s %20s" % ("ell/H", "H/ell", "gamma_A", "gamma_B", "kI_A", "kI_B"))
    rows = []
    for eoh in SAMPLES_ELL_OVER_H:
        HoL = 1.0 / eoh
        gA, gB = gamma_A_analytic(HoL), gamma_B_analytic(HoL)
        kA, kB = kI_from_gamma(gA), kI_from_gamma(gB)
        # numerical routes
        solH = strip_bvp_collocation(HoL / 2.0)
        sol2H = strip_bvp_collocation(HoL)  # a = (2H)/(2 ell)
        gA_num, _ = gamma_num_from_solution(solH, HoL)
        gB_num = 2.0 * gA_num - gamma_num_from_solution(sol2H, 2.0 * HoL)[0]
        # FD route (fine grid: tridiagonal solve is cheap; error ~ O(N^-2))
        NFD = 40000
        etaH, vH = strip_bvp_fd(HoL / 2.0, NFD)
        eta2H, v2H = strip_bvp_fd(HoL, NFD)
        gA_fd = -integrate_energy(etaH, vH, np.gradient(vH, etaH, edge_order=2))
        gB_fd = 2.0 * gA_fd + integrate_energy(eta2H, v2H, np.gradient(v2H, eta2H, edge_order=2))
        say("%9.4f %9.6f %20.16f %20.16f %20.16f %20.16f" % (eoh, HoL, gA, gB, kA, kB))
        rows.append(
            dict(
                ell_over_H=eoh,
                H_over_ell=HoL,
                gamma_A=gA,
                gamma_B=gB,
                kI_A=kA,
                kI_B=kB,
                gamma_A_numerical_bvp=gA_num,
                gamma_B_numerical_bvp=gB_num,
                gamma_A_numerical_fd=gA_fd,
                gamma_B_numerical_fd=gB_fd,
            )
        )
    say("")
    tol_prod = 5e-9
    worst = 0.0
    for r in rows:
        for key, ref in (("gamma_A_numerical_bvp", r["gamma_A"]), ("gamma_B_numerical_bvp", r["gamma_B"]),
                         ("gamma_A_numerical_fd", r["gamma_A"]), ("gamma_B_numerical_fd", r["gamma_B"])):
            worst = max(worst, abs(r[key] - ref))
    assert worst < tol_prod, "numerical vs analytic mismatch %g exceeds %g" % (worst, tol_prod)
    say("[MAIN] max |numerical - analytic| over all samples and routes: %.3e < %g  -> PASS" % (worst, tol_prod))
    say("")

    # ---- controls (deliberately wrong) ---------------------------------------
    say("=" * 100)
    say("[CONTROLS] deliberately WRONG variants; actual outputs preserved; NOT production answers")
    ctrl = {}
    # C1: omit foundation energy (k/2) u^2  -> integrand v' + v'^2/2 only
    c1_vals = []
    for eoh in SAMPLES_ELL_OVER_H:
        HoL = 1.0 / eoh
        a = HoL / 2.0
        sol = strip_bvp_collocation(a)
        eta = np.linspace(0.0, a, 20001)
        v = sol.sol(eta)[0]
        vp = sol.sol(eta, 1)[0]
        I_c1 = -2.0 * v[0] + np.sum(simpson_weights(eta) * vp**2)  # v' term + elastic only
        c1_vals.append(-I_c1)  # wrong gamma_A
    ctrl["C1_omit_foundation_energy_gamma_A_wrong"] = c1_vals
    say("  C1 (omit foundation energy): wrong gamma_A at ell/H=%s -> %s"
        % (SAMPLES_ELL_OVER_H, ["%.6f" % c for c in c1_vals]))
    say("     correct gamma_A = %s ; max deviation %.3e  -> CONTROL FAILS AS EXPECTED"
        % (["%.6f" % gamma_A_analytic(1.0 / e) for e in SAMPLES_ELL_OVER_H],
           max(abs(c - gamma_A_analytic(1.0 / e)) for c, e in zip(c1_vals, SAMPLES_ELL_OVER_H))))
    # C2: incremental traction on faces instead of total traction: v'(0) = 0 -> v = 0
    a = 0.5  # S/ell = 1 (ell/H = 1 row)
    sol0 = solve_bvp(
        lambda x, y: np.vstack((y[1], y[0])),
        lambda ya, yb: np.array([ya[1], yb[0]]),  # v'(0) = 0 (WRONG: incremental traction)
        np.linspace(0.0, a, 32),
        np.zeros((2, 32)),
        tol=1e-12,
    )
    gA_c2, I_c2 = gamma_num_from_solution(sol0, 2 * a)
    ctrl["C2_incremental_traction_gamma_A_wrong"] = gA_c2
    say("  C2 (incremental instead of total traction on faces): collapses to v = 0,")
    say("     wrong gamma_A = %.3e (correct %.6f); max|v| = %.3e  -> CONTROL FAILS AS EXPECTED"
        % (gA_c2, gamma_A_analytic(1.0), float(np.max(np.abs(sol0.sol(np.linspace(0, a, 100))[0])))))
    # C3: swap initial/final strip spacings in Case B (ahead = H-wake, behind = 2H-wake)
    c3_vals = []
    for eoh in SAMPLES_ELL_OVER_H:
        HoL = 1.0 / eoh
        gA_num_H, _ = gamma_num_from_solution(strip_bvp_collocation(HoL / 2.0), HoL)
        gA_num_2H, _ = gamma_num_from_solution(strip_bvp_collocation(HoL), 2.0 * HoL)
        c3_vals.append(2.0 * gA_num_2H - gA_num_H)  # swapped states: 2*G'(2H) - G'(H)
    ctrl["C3_swapped_spacings_gamma_B_wrong"] = c3_vals
    say("  C3 (initial/final spacings swapped: 2*G'(2H) - G'(H)): wrong gamma_B at ell/H=%s -> %s"
        % (SAMPLES_ELL_OVER_H, ["%.6f" % c for c in c3_vals]))
    say("     correct gamma_B = %s  -> CONTROL FAILS AS EXPECTED"
        % (["%.6f" % gamma_B_analytic(1.0 / e) for e in SAMPLES_ELL_OVER_H]))
    # assertions that controls fail by a clear margin
    assert all(abs(c - gamma_A_analytic(1.0 / e)) > 1e-3 for c, e in zip(c1_vals, SAMPLES_ELL_OVER_H))
    assert abs(gA_c2) < 1e-9 and abs(gA_c2 - gamma_A_analytic(1.0)) > 1e-3
    assert all(abs(c - gamma_B_analytic(1.0 / e)) > 1e-3 for c, e in zip(c3_vals, SAMPLES_ELL_OVER_H)), \
        "C3 control unexpectedly close to correct"
    say("")

    # ---- write reference.json ------------------------------------------------
    doc = {
        "meta": {
            "task": "Energy release for periodic cracks in an elastic film (PROBLEM.md)",
            "source": "Z.C. Xia & J.W. Hutchinson, Crack patterns in thin films, JMPS 48 (2000) 1107-1131; Sec. 4.2; Eqs. 28, 30 (see also Eqs. 4-12, 27); Eq. 6 notation ambiguity noted in REFERENCE.md Sec. 7",
            "nu": NU,
            "fixed_variables": "h, E, nu, k, sigma0 fixed during each increment of crack advance (PROBLEM.md); nu = 0.3 for the tabulated values; ell/H values are normalized sampling locations, not material data",
            "definitions": {
                "ell": "sqrt(Ebar*h/k) wake decay (interaction) length; Ebar = E/(1-nu^2)",
                "gamma_A": "G_A*Ebar/(sigma0^2*ell) = tanh(H/(2*ell))  [dimensionless]",
                "gamma_B": "G_B*Ebar/(sigma0^2*ell) = 2*tanh(H/(2*ell)) - tanh(H/ell) = 2*G'(H) - G'(2H)  [dimensionless]",
                "kI_A": "KI_A/(sigma0*sqrt(ell)) = sqrt((1-nu^2)*gamma_A)  [dimensionless]",
                "kI_B": "KI_B/(sigma0*sqrt(ell)) = sqrt((1-nu^2)*gamma_B)  [dimensionless]",
            },
            "units": {
                "G": "J/m^2 (energy per advancing crack area; release per unit in-plane advance per tip divided by h)",
                "KI": "Pa*sqrt(m)",
                "sigma_c": "sigma_c = sqrt(Gc*Ebar/(ell*gamma(H))) [Pa]",
                "ell": "m; k in Pa/m",
            },
            "assumptions": [
                "plane-stress film on Winkler foundation f = -k u; linear elasticity; quasistatic",
                "u1 = 0 far from tips (mode-I symmetry); 1D wake states exactly solve the model there",
                "total traction sigma.n = 0 on crack faces (Delta_sigma.n = -sigma0 n), NOT incremental traction",
                "steady simultaneous advance (Case A, one moving tip per period H) vs steady sequential midway advance (Case B, one moving new tip per period 2H; old cracks fixed)",
                "pure mode I by reflection symmetry: KII = 0; G = KI^2/E",
                "energy from the problem's stated incremental functional h*sigma0*eps_gg + (h/2)*Delta_sigma*eps + (k/2)u^2",
            ],
            "limits": {
                "H/ell->0": "gamma_A ~ (H/2*ell)(1-(H/ell)^2/12), gamma_B ~ (H/ell)^3/4; sigma_c,A ~ sqrt(2*Gc*Ebar/H), sigma_c,B ~ 2*sqrt(Gc*Ebar)*ell/H^(3/2) (shielding; sequential cracking arrests at spacing of order ell)",
                "H/ell->inf": "gamma_A, gamma_B -> 1; sigma_c -> sqrt(Gc*Ebar/ell) (isolated semi-infinite crack, XH2000 Eqs. 10-12, 27); interaction negligible beyond H ~ 3*ell",
            },
            "numerical_validation": {
                "routes": "analytic; scipy.solve_bvp collocation strip BVP; independent 2nd-order finite-difference BVP; Simpson quadrature of the stated functional using only numerical v",
                "max_bvp_vs_exact": max_bvp_err,
                "max_energy_residual": max_I_err,
                "max_prod_mismatch": worst,
                "recommended_comparison_tolerances": {
                    "gamma": "1e-6 absolute (observed route agreement <= %.1e)" % worst,
                    "kI": "1e-6 absolute (sqrt roughly halves relative error)",
                    "G_dimensional": "1e-6 relative, plus any physical-input uncertainty (substrate wavelength dependence, k calibration) which dominates in practice",
                },
            },
            "controls_deliberately_wrong": ctrl,
            "grade": "none (reference preparation only; no experimental grading or physical validation implied)",
        },
        "rows": rows,
    }
    out = os.path.join(HERE, "reference.json")
    with open(out, "w") as f:
        json.dump(doc, f, indent=2)
    say("[OUT] wrote %s" % out)
    say("[DONE] all production checks passed; all three controls failed as expected.")
    with open(os.path.join(HERE, "reproduce_output.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
