import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, BilayerCrack, solve, crack_opening, interface_tractions
from central_crack import end_loads_uniform

nu = 1/3
Db = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]
Dh = [plane_strain_D(1.0, nu), plane_strain_D(1.0, nu)]

# BilayerCrack (crack along all x<0) with uniform tension on LEFT end only.
mesh = BilayerCrack(1.0, 2.0, 30.0, 3e-4, 2e-4)
f = end_loads_uniform(mesh, 1.0)
print("left-end total Fx:", f[0::2].sum())   # expect -(h+H) = -3
u = solve(mesh, Db, f)
rd, dd = crack_opening(mesh, u)
m = np.argmin(np.abs(rd-1.0))
print(f"bilayer crack w/ uniform left load: delta(x=-1) = {dd[m]:.4f}")
T = interface_tractions(mesh, Db, u, np.array([0.02]))
print("T(0.02):", T)

# Now check equilibrium of that state: reaction at right end? solve() pins only
# 3 dofs; unbalanced load -> huge drift. Print max|u|:
print("max|u|:", np.abs(u).max(), "(unbalanced load -> large rigid drift ok)")
