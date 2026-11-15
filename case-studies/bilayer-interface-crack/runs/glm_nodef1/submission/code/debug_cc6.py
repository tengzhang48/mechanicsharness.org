import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D
from central_crack import CentralCrack, end_loads_uniform
from scipy.sparse.linalg import spsolve

# test 1: homogeneous material
nu = 1/3
Dh = [plane_strain_D(1.0, nu)]*2
mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
f = end_loads_uniform(mesh, 1.0)
K = fem.assemble(mesh.coords, mesh.elems, mesh.mats, Dh)
pin = [2*mesh.nid[0, 0], 2*mesh.nid[0, 0]+1, 2*mesh.nid[0, mesh.ny-1]]
free = np.setdiff1d(np.arange(2*len(mesh.coords)), pin)
u = np.zeros(2*len(mesh.coords))
u[free] = spsolve(K[free][:, free], f[free])
rd, dd = fem.crack_opening(mesh, u)
m = np.argmin(np.abs(rd-1.0))
print(f"homogeneous: max|u|={np.abs(u).max():.3f} delta(x=-1)={dd[m]:.4e}"
      f" (expected ~ 2.2)")
rs = np.array([0.01, 0.05])
T = fem.interface_tractions(mesh, Dh, u, rs)
print("  T:", T, " K-ish:", T*np.sqrt(2*np.pi*rs), " (sigma sqrt(pi a)=1.7725)")

# test 2: bimaterial, displacement control (prescribe ux=d on right end)
Db = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]
mesh2 = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
K2 = fem.assemble(mesh2.coords, mesh2.elems, mesh2.mats, Db)
right = []
for j in range(mesh2.ny):
    right.append(mesh2.nid[mesh2.nx-1, j])
    if j < mesh2.ny-1:
        right.append(mesh2.vnid[mesh2.nx-1, j])
right = np.array(sorted(right))
A = mesh2.nid[0, 0]
pin = np.concatenate([[2*A, 2*A+1], 2*right])       # ux,uy@A, ux=d on right end
d0 = 1.0
f2 = np.zeros(2*len(mesh2.coords))
u2 = np.zeros(2*len(mesh2.coords))
free = np.setdiff1d(np.arange(2*len(mesh2.coords)), pin)
# split pin into fixed-zero and fixed-d:
uzero = np.concatenate([[2*A, 2*A+1]])
ud = 2*right
u2[uzero] = 0.0
u2[ud] = d0
rhs = f2 - K2 @ u2*0  # solve with u prescribed via reduction:
free_all = np.setdiff1d(np.arange(2*len(mesh2.coords)), np.concatenate([uzero, ud]))
u2[free_all] = spsolve(K2[free_all][:, free_all], -K2[free_all][:, np.concatenate([uzero, ud])] @ u2[np.concatenate([uzero, ud])])
rd2, dd2 = fem.crack_opening(mesh2, u2)
m2 = np.argmin(np.abs(rd2-1.0))
print(f"bimaterial disp-control: max|u|={np.abs(u2).max():.3f} "
      f"delta(x=-1)={dd2[m2]:.4e} (expected ~ 0.05-0.2)")
