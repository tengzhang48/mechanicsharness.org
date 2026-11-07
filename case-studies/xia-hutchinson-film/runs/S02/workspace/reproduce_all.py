import numpy as np
# Driver that reproduces the reported table. Reuses the saved FEM modules.
import fem_energy2 as feA   # Case A energy-difference FEM
import femb2 as feB         # Case B half-domain FEM

print("=== Case A (fem_energy2.solve_tip; sh>0 shorter -> use -sh) ===")
print("ell/H  G_A       gamma_A   KI_A/(s0 sqrt ell)")
for r in [0.5, 1.0, 1.5, 2.0]:
    ell = r
    nx, ny, Xb, Xf, sh = 900, 120, 25.0, 10.0, 10
    P0, dx = feA.solve_tip(ell, nx, ny, Xb, Xf, 0)
    P1, _  = feA.solve_tip(ell, nx, ny, Xb, Xf, -sh)   # longer crack
    G = -(P1 - P0) / (sh * dx) / 1.0
    gam = G * (1.0 / (1 - 0.3 ** 2)) / ell
    print(f"{r:5.2f}  {G:.6f}  {gam:.6f}  {np.sqrt(G / ell):.6f}")

print("\n=== Case B (femb2.solve_tip; half domain, factor 2 for full 2H) ===")
print("ell/H  G_B       gamma_B   KI_B/(s0 sqrt ell)")
for r in [0.5, 1.0, 1.5, 2.0]:
    ell = r
    nx, ny, Xb, Xf, sh = 1000, 90, 25.0, 10.0, 8
    P0, dx = feB.solve_tip(ell, nx, ny, Xb, Xf, 0)
    P1, _  = feB.solve_tip(ell, nx, ny, Xb, Xf, sh)    # longer crack
    G = 2.0 * (-(P1 - P0) / (sh * dx) / 1.0)
    gam = G * (1.0 / (1 - 0.3 ** 2)) / ell
    print(f"{r:5.2f}  {G:.6f}  {gam:.6f}  {np.sqrt(G / ell):.6f}")