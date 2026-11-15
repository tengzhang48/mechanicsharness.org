import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, solve, crack_opening
from central_crack import CentralCrack, end_loads_uniform
from scipy.sparse.linalg import splu

nu = 1/3
Db = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]
mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
f = end_loads_uniform(mesh, 1.0, both=True)
K = fem.assemble(mesh.coords, mesh.elems, mesh.mats, Db).tocsc()

# pin ux,uy at bottom-left and ux,uy at top-left AND uy at right-bottom
A = mesh.nid[0, 0]; B = mesh.nid[0, mesh.ny-1]; C = mesh.nid[mesh.nx-1, 0]
pin = np.array([2*A, 2*A+1, 2*B, 2*B+1, 2*C+1])
free = np.setdiff1d(np.arange(2*len(mesh.coords)), pin)
u = np.zeros(2*len(mesh.coords))
u[free] = splu(K[free][:, free].tocsc()).solve(f[free])
rd, dd = crack_opening(mesh, u)
m = np.argmin(np.abs(rd-1.0))
print(f"heavy pinning: max|u|={np.abs(u).max():.4f} delta(x=-1)={dd[m]:.5f}")
i = np.searchsorted(mesh.X, 20.0)
n = mesh.nid[i, mesh.j0]
print("u at (20, 0):", u[2*n], u[2*n+1])
