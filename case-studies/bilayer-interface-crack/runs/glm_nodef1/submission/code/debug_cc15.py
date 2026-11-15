import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, solve
from central_crack import CentralCrack, end_loads_uniform

nu = 1/3
Db = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]
mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
print("X span:", mesh.X[0], mesh.X[-1], "nx:", mesh.nx)
f = end_loads_uniform(mesh, 1.0, both=True)
u = solve(mesh, Db, f)
K = fem.assemble(mesh.coords, mesh.elems, mesh.mats, Db)
print("||Ku-f||/||f||:", np.linalg.norm(K@u-f)/np.linalg.norm(f))
pair = min(mesh.pairs, key=lambda p: abs(p[2]+1.0))
b, d, x1 = pair
print("pair ids:", b, d, " coords:", mesh.coords[b], mesh.coords[d])
print("u_upper:", u[2*b], u[2*b+1])
print("u_lower:", u[2*d], u[2*d+1])
# are b and d distinct dof entries in K?
print("K row nnz upper:", K[2*b].nnz, " lower:", K[2*d].nnz)
print("K[2b,2d]:", K[2*b, 2*d], " K[2b+1,2d+1]:", K[2*b+1, 2*d+1])
# which elements use d?
ed = np.where((mesh.elems == d).any(axis=1))[0]
print("elements containing lower dup d:", ed)
for e in ed[:4]:
    print("   elem", e, "node coords:", np.round(mesh.coords[mesh.elems[e]], 4).tolist())
