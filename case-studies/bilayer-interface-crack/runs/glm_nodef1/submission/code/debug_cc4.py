import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D
from central_crack import CentralCrack, end_loads_uniform
from scipy.sparse.linalg import spsolve

nu = 1/3
D = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]
mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
f = end_loads_uniform(mesh, 1.0)
K = fem.assemble(mesh.coords, mesh.elems, mesh.mats, D)
A = mesh.nid[0, 0]; B = mesh.nid[0, mesh.ny-1]
C = mesh.nid[mesh.nx-1, 0]
for pin, label in [
    ([2*A, 2*A+1, 2*B], "ux,uy@A + ux@B  (current fem.solve)"),
    ([2*A, 2*A+1, 2*B+1], "ux,uy@A + uy@B"),
    ([2*A, 2*A+1, 2*C+1], "ux,uy@A + uy@C(right-bottom)"),
]:
    free = np.setdiff1d(np.arange(2*len(mesh.coords)), pin)
    u = np.zeros(2*len(mesh.coords))
    u[free] = spsolve(K[free][:, free], f[free])
    r = np.linalg.norm(K @ u - f)/np.linalg.norm(f)
    # crack opening at x=-1
    rd, dd = fem.crack_opening(mesh, u)
    m = np.argmin(np.abs(rd-1.0))
    print(f"{label}: res={r:.1e} max|u|={np.abs(u).max():.3f} "
          f"delta(x=-1)={dd[m]:.3e}")
