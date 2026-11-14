import numpy as np, sys
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem

p = WavyFrontProblem(A=0.0, harm=1, Lcr=2.0, Lx=3.5, Nxf=2048, Nz=32,
                     M=10, Nzmode=2, taper=(0.30,0.95))
W, u, sig = p._window_u()
print("W stats: max", W.max(), "frac nonzero", (W>0).mean())
print("u stats:", u.min(), u.max())
g = p.basis_grid(1, 0, 0, W, u, sig)
print("basis m0 |g| max:", np.abs(g[1]).max())
ii, jj = p.collocation(margin_dx=15.0)
print("colloc pts:", len(ii), " s range:", (p.a[0,jj]-p.X[ii,jj]).min(),
      (p.a[0,jj]-p.X[ii,jj]).max())
info = p.solve(ridge=0.0)
print("solve:", info)
print("c norm:", np.linalg.norm(p.c))
th, K1, K2, K3 = p.K_along_front()
print("K1:", K1.mean())
# check matrix conditioning quickly
