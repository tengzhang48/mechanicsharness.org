import numpy as np
import sys, warnings
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, solve, crack_opening, interface_tractions
from central_crack import CentralCrack, end_loads_uniform

nu = 1/3
Db = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]

# variant: a = 20 -> crack duplicated for all x<0 (like BilayerCrack)
mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 20.0)
print("pairs:", len(mesh.pairs), " X range:", mesh.X[0], mesh.X[-1])
sigma = {1: lambda x2: -0.0 - 12*1.0*(x2 - 0.5),
         2: lambda x2:  0.0/2 + 12*(1.0 + 0.0)*(x2 + 1.0)}
f = fem.end_loads(mesh, sigma)
with warnings.catch_warnings():
    warnings.simplefilter("error")
    try:
        u = solve(mesh, Db, f)
    except Exception as ex:
        print("solve raised:", type(ex).__name__, ex)
        warnings.simplefilter("ignore")
        u = solve(mesh, Db, f)
rd, dd = crack_opening(mesh, u)
m = (rd > 0.05) & (rd < 0.15)
print("delta near x=-0.1:", dd[m][:3])
T = interface_tractions(mesh, Db, u, np.array([0.02, 0.05]))
print("T(0.02):", T, " T*sqrt(2 pi r):", T*np.sqrt(2*np.pi*np.array([0.02, 0.05])))
print("(BilayerCrack case-B reference K ~ 1.62+0.12j)")
