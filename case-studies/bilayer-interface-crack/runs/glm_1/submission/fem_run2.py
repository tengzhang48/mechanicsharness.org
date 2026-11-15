import numpy as np
from fem import solve_fem

print("TARGET, mesh refinement + fit-range study:")
for tag, P, M in [("A", 1.0, 0.0), ("B", 0.0, 1.0), ("AB", 1.0, 1.0)]:
    for dmin in [2.5e-4, 1e-4, 4e-5]:
        row = f"  case {tag} dmin={dmin:g}: "
        for (rmin, rmax) in [(8*dmin, 0.15), (8*dmin, 0.25), (4*dmin, 0.10)]:
            r = solve_fem(1.0, 2.0, 40.0, 9.0, 1.0, 1/3, P, M, dmin=dmin,
                          ny1=14, ny2=22, fit_rmin=rmin, fit_rmax=rmax)
            row += (f"[{rmin:.4g}-{rmax:.2g}] |K|={r['Kabs']:.4f} psi={r['psi_deg']:+.3f}  ")
        print(row)
    print()
