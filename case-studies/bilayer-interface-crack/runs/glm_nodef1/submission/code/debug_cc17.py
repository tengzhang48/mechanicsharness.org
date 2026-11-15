import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, solve
from central_crack import CentralCrack, end_loads_uniform

nu = 1/3
Db = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]
mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
f = end_loads_uniform(mesh, 1.0, both=True)
u = solve(mesh, Db, f)

e = mesh.elem_of_cell[mesh.j0][np.searchsorted(mesh.X, 0.5)]
print("element near x=0.5, upper layer:")
for k in range(8):
    n = e[k]
    print(f"  node {n:6d} ({mesh.coords[n][0]:+8.4f},{mesh.coords[n][1]:+.4f})  "
          f"ux={u[2*n]:+12.5f} uy={u[2*n+1]:+12.5f}")
s = fem.elem_stress_at(mesh.coords, e, Db[0], u, 0.5, 0.0001)
print("stress at (0.5, 0.0001):", s)

e2 = mesh.elem_of_cell[mesh.j0-1][np.searchsorted(mesh.X, 0.5)]
s2 = fem.elem_stress_at(mesh.coords, e2, Db[1], u, 0.5, -0.0001)
print("lower stress at (0.5, -0.0001):", s2)

# look at vertical column at x=+10: how does uy vary with y?
i = np.searchsorted(mesh.X, 10.0)
print("\ncolumn at x~10 (all rows):")
for j in range(0, mesh.ny, 3):
    n = mesh.nid[i, j]
    print(f"  y={mesh.Y[j]:+8.4f}  ux={u[2*n]:+11.4f}  uy={u[2*n+1]:+11.4f}")
