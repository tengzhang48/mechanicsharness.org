import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, assemble
from central_crack import CentralCrack, end_loads_uniform
from scipy.sparse.linalg import spsolve

nu = 1/3
D = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]
mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
f = end_loads_uniform(mesh, 1.0)
K = assemble(mesh.coords, mesh.elems, mesh.mats, D)
pin = [(mesh.nid[0, 0], 'both'), (mesh.nid[0, mesh.ny-1], 'y')]
fixed = [2*mesh.nid[0, 0], 2*mesh.nid[0, 0]+1, 2*mesh.nid[0, mesh.ny-1]]
alldof = np.arange(2*len(mesh.coords))
free = np.setdiff1d(alldof, fixed)
u = np.zeros(2*len(mesh.coords))
u[free] = spsolve(K[free][:, free], f[free])
r = K @ u - f
print("residual relative:", np.linalg.norm(r)/np.linalg.norm(f))
print("reactions at pinned dofs:", r[fixed])
# check: any node not in any element?
inuse = np.zeros(len(mesh.coords), bool)
inuse[mesh.elems.ravel()] = True
print("nodes not in any element:", np.where(~inuse)[0])
