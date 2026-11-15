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
j0 = mesh.j0
print("upper interface row nodes (x, uy, ux):")
for i in range(mesh.nx):
    n = mesh.nid[i, j0]
    x = mesh.X[i]
    if -3.5 < x < 3.5 and (i % 7 == 0 or abs(abs(x)-2) < 0.01 or abs(x) < 0.001):
        print(f"   x={x:+8.4f}  ux={u[2*n]:+10.4f}  uy={u[2*n+1]:+10.4f}")
print()
print("row at y=-H (bottom surface), every ~10th node:")
for i in range(0, mesh.nx, 12):
    n = mesh.nid[i, 0]
    print(f"   x={mesh.X[i]:+8.3f}  ux={u[2*n]:+10.4f}  uy={u[2*n+1]:+10.4f}")
