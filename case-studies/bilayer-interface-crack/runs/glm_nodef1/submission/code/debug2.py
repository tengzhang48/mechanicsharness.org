import numpy as np
import sys
sys.path.insert(0, '/work/code')
from fem import (BilayerCrack, plane_strain_D, end_loads, solve,
                 interface_tractions, crack_opening)
from fits import fit_series

nu = 1/3
D = [plane_strain_D(1.0, nu)]*2
h = H = 1.0
L = 30.0
P, M = 0.0, 1.0
sigma = {1: lambda x2: -P/h - 12*M*(x2 - h/2)/h**3,
         2: lambda x2:  P/H + 12*(M + P*(h+H)/2)*(x2 + H/2)/H**3}

rs = np.geomspace(0.002, 0.2, 25)
print("pointwise T*sqrt(2 pi r) at fixed r for different meshes:")
print("r        " + "  ".join(f"d1={d1:.0e}/{d2:.0e}" for d1, d2 in
      [(3e-4, 5e-3), (1.5e-4, 5e-3), (3e-4, 2.5e-3), (7.5e-5, 1.25e-3), (1.5e-4, 1.25e-3)]))
data = {}
for (d1, d2) in [(3e-4, 5e-3), (1.5e-4, 5e-3), (3e-4, 2.5e-3), (7.5e-5, 1.25e-3), (1.5e-4, 1.25e-3)]:
    mesh = BilayerCrack(h, H, L, d1, d2)
    f = end_loads(mesh, sigma)
    u = solve(mesh, D, f)
    data[(d1, d2)] = interface_tractions(mesh, D, u, rs)*np.sqrt(2*np.pi*rs)
for k, r in enumerate(rs):
    row = f"{r:.4f}   "
    for key in data:
        row += f"  {data[key][k].real:9.4f}"
    print(row)
print("\nexact K_I = 3.464102 (pure real)")
print("\nfits K (basis 1, r^0.5, r, r^1.5) over sub-ranges, mesh d1=1.5e-4/d2=1.25e-3:")
mesh = BilayerCrack(h, H, L, 1.5e-4, 1.25e-3)
f = end_loads(mesh, sigma)
u = solve(mesh, D, f)
Tn = interface_tractions(mesh, D, u, rs)*np.sqrt(2*np.pi*rs)
for rmax in [0.02, 0.05, 0.1, 0.2]:
    row = f"rmax={rmax:5.2f}: "
    for rmin in [0.002, 0.005, 0.01]:
        m = (rs >= rmin) & (rs <= rmax)
        c = fit_series(rs[m], Tn[m], (0.0, 0.5, 1.0, 1.5))
        row += f"  rmin={rmin:6.3f}: K={c[0].real:9.4f}{c[0].imag:+.4f}j"
    print(row)
print("exact K_I = 3.464102")
