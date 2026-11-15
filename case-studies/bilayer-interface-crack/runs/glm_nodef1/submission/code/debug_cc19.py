import numpy as np
import sys, warnings
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, crack_opening
from central2 import Central2, X, a, L, D

mesh = Central2(1.0, 2.0, L, 4e-4, 2e-4, a, X)
f = np.zeros(2*len(mesh.coords))
for edges, sgn in ((mesh.left_edges, -1.0), (mesh.right_edges, +1.0)):
    for (n, yb, yt, mat) in edges:
        half = 0.5*(yt-yb)
        for (t, w) in ((-np.sqrt(3/5), 5/9), (0.0, 8/9), (np.sqrt(3/5), 5/9)):
            Na = 0.5*t*(t-1); Nb = 0.5*t*(t+1); Nm = 1-t*t
            fx = sgn*1.0*half*w
            f[2*n[0]] += Na*fx
            f[2*n[1]] += Nb*fx
            f[2*n[2]] += Nm*fx
K = fem.assemble(mesh.coords, mesh.elems, mesh.mats, D).tocsc()
Apin = mesh.nid[0, 0]
pin = [2*Apin, 2*Apin+1, 2*mesh.nid[0, mesh.ny-1]]
free = np.setdiff1d(np.arange(2*len(mesh.coords)), pin)
Kf = K[free][:, free].tocsc()
ff = f[free]
print("cond estimate via diagonal scaling:")
d = Kf.diagonal()
print("  min/max diag:", d.min(), d.max())
from scipy.sparse.linalg import splu, spsolve, cg, LinearOperator
for name, method in [("splu-COLAMD", lambda: splu(Kf).solve(ff)),
                     ("splu-NATURAL", lambda: splu(Kf, permc_spec="NATURAL").solve(ff)),
                     ("splu-MMD_AT_PLUS_A", lambda: splu(Kf, permc_spec="MMD_AT_PLUS_A").solve(ff)),
                     ("cg", lambda: cg(Kf, ff, rtol=1e-10, maxiter=20000)[0])]:
    try:
        x = method()
        u = np.zeros(2*len(mesh.coords)); u[free] = x
        rd, dd = crack_opening(mesh, u)
        mm = np.argmin(np.abs(rd-1.0))
        print(f"  {name}: max|u|={np.abs(u).max():10.3f} delta={dd[mm]:+.4e} "
              f"resid={np.linalg.norm(Kf@x-ff):.2e}")
    except Exception as e:
        print(f"  {name}: FAILED {e}")
