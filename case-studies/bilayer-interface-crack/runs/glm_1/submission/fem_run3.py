import numpy as np
from fem import solve_fem

def run_all(P, M, tag):
    out = []
    for (dmin, L, ny1, ny2, rmax, nt) in [
        (1e-4, 40.0, 14, 22, 0.15, 3),
        (1e-4, 40.0, 14, 22, 0.20, 4),
        (1e-4, 60.0, 18, 28, 0.15, 3),
        (4e-5, 40.0, 14, 22, 0.15, 3),
        (4e-5, 40.0, 16, 26, 0.15, 4),
        (4e-5, 60.0, 18, 28, 0.15, 3),
        (4e-5, 60.0, 18, 28, 0.10, 3),
        (1.5e-5, 60.0, 18, 28, 0.15, 3),
    ]:
        r = solve_fem(1.0, 2.0, L, 9.0, 1.0, 1/3, P, M, dmin=dmin,
                      ny1=ny1, ny2=ny2, fit_rmin=8 * dmin, fit_rmax=rmax,
                      nterms=nt)
        out.append(r['psi_deg'])
        print(f"  {tag} dmin={dmin:g} L={L:.0f} ny=({ny1},{ny2}) rmax={rmax} nt={nt}: "
              f"|K|={r['Kabs']:.4f} psi={r['psi_deg']:+.4f} ndof={r['ndof']}")
    return np.array(out)

print("case A:")
pA = run_all(1.0, 0.0, "A")
print("case B:")
pB = run_all(0.0, 1.0, "B")
print("case AB (P=1, M=1):")
pC = run_all(1.0, 1.0, "AB")
print()
print(f"medians: psi_A={np.median(pA):+.4f} (spread {pA.max()-pA.min():.3f})  "
      f"psi_B={np.median(pB):+.4f} (spread {pB.max()-pB.min():.3f})  "
      f"psi_AB={np.median(pC):+.4f} (spread {pC.max()-pC.min():.3f})")
print(f"psi_A - psi_B = {np.median(pA)-np.median(pB):+.3f}  (exact |.| = 47.9705)")
# combined prediction from exact magnitudes + FEM phases
KA, KB = 1.9331267849436742, 1.629800601300662
z = KA*np.exp(1j*np.radians(np.median(pA))) + KB*np.exp(1j*np.radians(np.median(pB)))
print(f"psi_AB from exact |K| + FEM phases = {np.degrees(np.angle(z)):+.4f}")
