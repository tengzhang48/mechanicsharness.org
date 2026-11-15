import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, solve, crack_opening
from central_crack import CentralCrack, end_loads_uniform

nu = 1/3
Db = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]
mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
f = end_loads_uniform(mesh, 1.0, both=True)
Fx = f[0::2]; y = mesh.coords[:, 1]
print("sum Fx:", Fx.sum(), "  sum(-y Fx) =", -(y*Fx).sum())
u = solve(mesh, Db, f)
print("max|u|:", np.abs(u).max())
i_mid1 = mesh.nid[np.searchsorted(mesh.X, 30.0), mesh.j0+3]
i_mid2 = mesh.nid[np.searchsorted(mesh.X, 30.0), mesh.j0-3]
print("u at (30, y~0.3):", u[2*i_mid1], u[2*i_mid1+1])
print("u at (30, y~-0.3):", u[2*i_mid2], u[2*i_mid2+1])
el = mesh.elem_of_cell[mesh.j0+3][np.searchsorted(mesh.X, 30.0)]
s = fem.elem_stress_at(mesh.coords, el, Db[0], u, 30.0, 0.3)
print("stress at (30, 0.3):", s)
el2 = mesh.elem_of_cell[mesh.j0-4][np.searchsorted(mesh.X, 30.0)]
s2 = fem.elem_stress_at(mesh.coords, el2, Db[1], u, 30.0, -0.3)
print("stress at (30, -0.3):", s2)
