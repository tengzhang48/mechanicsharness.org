"""
Final driver: validations + final numbers for results.json.

Summary of method
-----------------
1. eps: exact closed form.
2. G (energy release rate): exact, from the steady-state argument.  The
   prescribed far-behind stress fields are exact x-independent elasticity
   solutions (each layer separately in equilibrium, traction-free outer
   faces, zero interface tractions), and the bonded region far ahead is
   stress free.  For steady advance of the semi-infinite crack under the
   dead loads (P, M), G = (strain energy per unit length far behind)
   - (far ahead = 0) + load work - load work = U_behind:
       dW/da = 2 U_behind (dead loads), dU/da = U_behind - 0
       =>  G = dW/da - dU/da = U_behind.
   Verified against the classical moment-loaded DCB formula
   G = 12 M^2/(E' t^3) and by the FEM.
3. |K|: exact from G via the interface relation
       G = (Lam1 + Lam2) |K|^2 / (16 cosh^2(pi eps)),
   Lam_k = (kappa_k + 1)/mu_k  (homogeneous limit reproduces G = K^2/E').
4. Phase angles psi: from an independent FEM solution of the actual
   boundary-value problem (fem.py): bilinear-quads, tip-graded mesh,
   crack-face duplication on x<0, shared nodes on x>0, beam tractions at
   x = -L, traction-free elsewhere; K from a 3-term fit of the interface
   tractions ahead of the tip with the oscillatory factor e^{i eps ln r}.
   Cross-checks: symmetric DCB limit (psi = 0 exactly, |K| = sqrt(12)),
   homogeneous strip case A (|K| = sqrt(7)), convergence in mesh and fit
   window, linearity in (P, M), and the exact energy relation
   cos(psi_A - psi_B) = (G(1,1) - G_A - G_B)/(2 sqrt(G_A G_B)).

Units: h = 1 (H = 2), E2 = 1 (E1 = 9), nu = 1/3 for both layers.
"""
import json
import numpy as np

from sie_solve import Problem            # exact eps, G, |K| (energy)
from fem import solve_fem

LOG = []


def log(s=""):
    print(s)
    LOG.append(str(s))


def main():
    # ======================= exact quantities ==============================
    pr = Problem(1.0, 2.0, 9.0, 1.0, 1 / 3)
    GA, GB = pr.G_exact(1, 0), pr.G_exact(0, 1)
    G11 = pr.G_exact(1, 1)
    KA, KB = pr.Kabs_exact(1, 0), pr.Kabs_exact(0, 1)
    cosDpsi = (G11 - GA - GB) / (2 * np.sqrt(GA * GB))
    Dpsi_ex = np.degrees(np.arccos(cosDpsi))
    log("=== exact (energy) quantities ===")
    log(f"eps                = {pr.eps:.12f}")
    log(f"G_A E2 h / P^2     = {GA:.10f}")
    log(f"G_B E2 h^3 / M^2   = {GB:.10f}")
    log(f"|K|_A sqrt(h)/P    = {KA:.10f}")
    log(f"|K|_B h^3/2/M      = {KB:.10f}")
    log(f"energy cross term: |psi_A - psi_B| = {Dpsi_ex:.4f} deg (exact)")

    # ======================= validations ===================================
    log("")
    log("=== validation V1: symmetric DCB (E1=E2, h=H, M only) ===")
    log("    by symmetry psi = 0 exactly; |K| = sqrt(12) = 3.464102 (energy)")
    for dmin in [1e-4, 1.5e-5]:
        r = solve_fem(1.0, 1.0, 60.0, 1.0, 1.0, 1 / 3, 0.0, 1.0, dmin=dmin,
                      ny1=18, ny2=18, fit_rmin=8 * dmin, fit_rmax=0.15)
        log(f"    dmin={dmin:g}: |K|={r['Kabs']:.5f} "
            f"(err {abs(r['Kabs']/np.sqrt(12)-1)*100:+.2f}%), psi={r['psi_deg']:+.4f} deg")
    log("=== validation V2: homogeneous strip, load case A (mixed mode) ===")
    log("    |K| = sqrt(7) = 2.645751 (energy)")
    r = solve_fem(1.0, 1.0, 60.0, 1.0, 1.0, 1 / 3, 1.0, 0.0, dmin=1.5e-5,
                  ny1=18, ny2=18, fit_rmin=1.2e-4, fit_rmax=0.15)
    log(f"    |K|={r['Kabs']:.5f} (err {abs(r['Kabs']/np.sqrt(7)-1)*100:+.2f}%)")

    # ======================= target runs ===================================
    log("")
    log("=== target problem: FEM convergence (phases) ===")
    res = {}
    for tag, P, M in [("A", 1.0, 0.0), ("B", 0.0, 1.0), ("AB", 1.0, 1.0)]:
        psis, kabs = [], []
        for dmin in [1e-4, 4e-5, 1.5e-5]:
            r = solve_fem(1.0, 2.0, 60.0, 9.0, 1.0, 1 / 3, P, M, dmin=dmin,
                          ny1=18, ny2=28, fit_rmin=8 * dmin, fit_rmax=0.15)
            psis.append(r["psi_deg"])
            kabs.append(r["Kabs"])
            log(f"    case {tag:2s} dmin={dmin:g}: |K|={r['Kabs']:.5f} "
                f"psi={r['psi_deg']:+.4f}")
        res[tag] = (np.median(psis), np.array(kabs))
        log(f"    case {tag:2s}: psi median {np.median(psis):+.4f} "
            f"(spread {max(psis)-min(psis):.3f} deg)")

    psiA, psiB, psiC = res["A"][0], res["B"][0], res["AB"][0]
    log("")
    log("=== cross-checks ===")
    log(f"    FEM |K|_A vs exact: {res['A'][1][-1]:.5f} vs {KA:.5f} "
        f"({abs(res['A'][1][-1]/KA-1)*100:+.2f}%)")
    log(f"    FEM |K|_B vs exact: {res['B'][1][-1]:.5f} vs {KB:.5f} "
        f"({abs(res['B'][1][-1]/KB-1)*100:+.2f}%)")
    log(f"    FEM |K|_AB vs exact {np.sqrt(16*np.cosh(np.pi*pr.eps)**2*G11/pr.Lam_sum()):.5f}: "
        f"{res['AB'][1][-1]:.5f} ({abs(res['AB'][1][-1]/np.sqrt(16*np.cosh(np.pi*pr.eps)**2*G11/pr.Lam_sum())-1)*100:+.2f}%)")
    log(f"    psi_A - psi_B (FEM) = {psiA-psiB:+.3f} deg vs exact {Dpsi_ex:.4f} "
        f"(diff {abs(psiA-psiB-Dpsi_ex):.3f} deg)")
    z = KA * np.exp(1j * np.radians(psiA)) + KB * np.exp(1j * np.radians(psiB))
    log(f"    psi_both: direct FEM {psiC:+.3f} deg; exact-|K| + FEM phases "
        f"{np.degrees(np.angle(z)):+.3f} deg")
    # linearity of the FEM solution
    rA = solve_fem(1.0, 2.0, 40.0, 9.0, 1.0, 1 / 3, 1.0, 0.0, dmin=4e-5,
                   ny1=14, ny2=22, fit_rmin=3.2e-4, fit_rmax=0.15)
    rB = solve_fem(1.0, 2.0, 40.0, 9.0, 1.0, 1 / 3, 0.0, 1.0, dmin=4e-5,
                   ny1=14, ny2=22, fit_rmin=3.2e-4, fit_rmax=0.15)
    rAB = solve_fem(1.0, 2.0, 40.0, 9.0, 1.0, 1 / 3, 1.0, 1.0, dmin=4e-5,
                    ny1=14, ny2=22, fit_rmin=3.2e-4, fit_rmax=0.15)
    log(f"    linearity: |K(1,1)| = {rAB['Kabs']:.5f} vs |K_A+K_B| = "
        f"{abs(rA['K']+rB['K']):.5f}")

    # ======================= results =======================================
    out = {
        "eps": float(pr.eps),
        "caseA": {
            "G_E2_h_over_P2": float(GA),
            "Kabs_sqrt_h_over_P": float(KA),
            "psi_deg": float(round(psiA, 2)),
        },
        "caseB": {
            "G_E2_h3_over_M2": float(GB),
            "Kabs_h32_over_M": float(KB),
            "psi_deg": float(round(psiB, 2)),
        },
        "psi_both_deg": float(round(psiC, 2)),
    }
    log("")
    log(json.dumps(out, indent=2))
    with open("results.json", "w") as f:
        json.dump(out, f, indent=2)
    with open("run_log.txt", "w") as f:
        f.write("\n".join(LOG))
    log("wrote results.json, run_log.txt")


if __name__ == "__main__":
    main()
