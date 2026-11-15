import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D
from central_crack import CentralCrack
from scipy.sparse.linalg import spsolve

nu = 1/3
mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
D = [plane_strain_D(1.0, nu)]*2
K = fem.assemble(mesh.coords, mesh.elems, mesh.mats, D)
# take an upper crack-face node at x ~ -1 and its lower dup
pair = min(mesh.pairs, key=lambda p: abs(p[2]+1.0))
b, d, x1 = pair
print("probe pair:", b, d, x1, mesh.coords[b], mesh.coords[d])
# apply unit force in y at upper node; fix far nodes; see displacement of lower dup
A = mesh.nid[0, 0]
pin = [2*A, 2*A+1, 2*mesh.nid[mesh.nx-1, 0]]
free = np.setdiff1d(np.arange(2*len(mesh.coords)), pin)
f = np.zeros(2*len(mesh.coords))
f[2*b+1] = 1.0
u = np.zeros(2*len(mesh.coords))
u[free] = spsolve(K[free][:, free], f[free])
print("u(upper node):", u[2*b], u[2*b+1])
print("u(lower dup):", u[2*d], u[2*d+1])
print("-> if lower dup moves like upper, faces are coupled")
# also check: which elements contain node b and d?
eb = np.where((mesh.elems == b).any(axis=1))[0]
ed = np.where((mesh.elems == d).any(axis=1))[0]
print("elements with upper node:", len(eb), "y-centroids:", mesh.coords[mesh.elems[eb][:, 0]][:, 1])
print("elements with lower dup:", len(ed), "y-centroids:", mesh.coords[mesh.elems[ed][:, 0]][:, 1])
