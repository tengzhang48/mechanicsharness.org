import numpy as np
import sys
sys.path.insert(0, '/work/code')
from fem import (BilayerCrack, plane_strain_D, end_loads, solve,
                 interface_tractions, crack_opening)

nu = 1/3
D = [plane_strain_D(1.0, nu)]*2
h = H = 1.0
L = 30.0

mesh = BilayerCrack(h, H, L, 1e-3, 8e-3)
print("nodes:", len(mesh.coords), "elems:", len(mesh.elems))
print("dup pairs:", len(mesh.pairs), "first few:", mesh.pairs[:3])

P, M = 0.0, 1.0
sigma = {1: lambda x2: -P/h - 12*M*(x2 - h/2)/h**3,
         2: lambda x2:  P/H + 12*(M + P*(h+H)/2)*(x2 + H/2)/H**3}
f = end_loads(mesh, sigma)
print("total Fx =", f[0::2].sum(), " (should be ~0)")
# moment about origin of applied loads (x*Fy - y*Fx), nodes at mesh coords
mom = np.sum(mesh.coords[:, 1]*(-f[0::2]))  # only Fx contributes (Fy=0): M = -y Fx
print("total moment about x1=0 of end loads =", mom, " (should be ~0)")

u = solve(mesh, D, f)
rd, dd = crack_opening(mesh, u)
sel = (rd > 0.01) & (rd < 0.5)
print("delta samples (r, |delta|):")
for k in np.where(sel)[0][:8]:
    print(f"   r={rd[k]:.4f}  delta={dd[k]:+.6e}")

rs = np.geomspace(0.003, 0.12, 9)
T = interface_tractions(mesh, D, u, rs)
print("T(r) samples:")
for r, t in zip(rs, T):
    print(f"   r={r:.4f}  T={t:.5f}   T*sqrt(2 pi r)={t*np.sqrt(2*np.pi*r):.5f} (exact K=3.4641)")
