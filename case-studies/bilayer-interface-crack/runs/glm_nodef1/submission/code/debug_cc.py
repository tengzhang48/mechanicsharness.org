import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, solve, interface_tractions, crack_opening
from central_crack import CentralCrack, end_loads_uniform

nu = 1/3
D = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]
mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
print("dup pairs:", len(mesh.pairs), "first/last:", mesh.pairs[:2], mesh.pairs[-2:])
f = end_loads_uniform(mesh, 1.0)
print("total Fx:", f[0::2].sum())
u = solve(mesh, D, f)
print("max |u|:", np.abs(u).max())
# displacement at far right mid-height (should be ~ sigma*L/Ebar ~ 40/10 ~ 4-ish)
i = mesh.nid[mesh.nx-1, mesh.j0]
print("u at right end interface:", u[2*i], u[2*i+1])
rd, dd = fem.crack_opening(mesh, u)
m = (rd > 0.5) & (rd < 1.5)
print("delta near crack middle (r, delta):")
for k in np.where(m)[0][:5]:
    print(f"   r={rd[k]:.3f} (x={-rd[k]:.3f})  delta={dd[k]:+.6e}")
rs = np.array([0.01, 0.05])
T = interface_tractions(mesh, D, u, rs)
print("T(r):", T)
# also check sigma11 near the right end (should be ~ 1)
el = mesh.elem_of_cell[mesh.j0][mesh.nx-3]
s = fem.elem_stress_at(mesh.coords, el, D[0], u, mesh.X[mesh.nx-4], 0.5)
print("far-field stress near right end:", s)
