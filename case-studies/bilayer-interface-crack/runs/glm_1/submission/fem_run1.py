import numpy as np
from fem import solve_fem

# --- validation 2: DCB strip, load case A: |K| = sqrt(7) (mixed mode, eps=0)
print("DCB case A (|K| exact sqrt(7) = 2.645751):")
for dmin in [1.5e-3, 6e-4, 2.5e-4, 1e-4]:
    r = solve_fem(1.0, 1.0, 30.0, 1.0, 1.0, 1/3, 1.0, 0.0, dmin=dmin)
    print(f"  dmin={dmin:g}: |K|={r['Kabs']:.5f} err={abs(r['Kabs']/np.sqrt(7)-1)*100:+.2f}% "
          f"psi={r['psi_deg']:+.4f} ndof={r['ndof']} fitresid={r['resid']:.1e}")

# --- target problem ---------------------------------------------------------
print("TARGET (E1=9, E2=1, nu=1/3, h=1, H=2):")
for tag, P, M in [("A", 1.0, 0.0), ("B", 0.0, 1.0)]:
    for dmin in [1.5e-3, 6e-4, 2.5e-4, 1e-4]:
        r = solve_fem(1.0, 2.0, 40.0, 9.0, 1.0, 1/3, P, M, dmin=dmin,
                      ny1=14, ny2=22)
        print(f"  case {tag} dmin={dmin:g}: |K|={r['Kabs']:.5f} psi={r['psi_deg']:+.4f} "
              f"ndof={r['ndof']} fitr={r['fit_r'][0]:.4f}..{r['fit_r'][1]:.2f} "
              f"fitresid={r['resid']:.1e}")
