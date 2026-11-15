import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, solve, crack_opening, interface_tractions, BilayerCrack
from central_crack import CentralCrack

nu = 1/3
Db = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]
h, H = 1.0, 2.0

sigma = {1: lambda x2: -12*(x2 - h/2),
         2: lambda x2: 12*(x2 + H/2)/H**3}

for tag, mesh in [("Bilayer", BilayerCrack(h, H, 10.0, 1e-3, 5e-4)),
                  ("Central a=5", CentralCrack(h, H, 10.0, 1e-3, 5e-4, 5.0))]:
    f = fem.end_loads(mesh, sigma)
    u = solve(mesh, Db, f)
    rd, dd = crack_opening(mesh, u)
    m = np.argmin(np.abs(rd-0.05))
    T = interface_tractions(mesh, Db, u, np.array([0.02]))
    print(f"{tag}: max|u|={np.abs(u).max():.4f} delta(r=0.05)={dd[m]:.4f} "
          f"T(0.02)*sqrt(2pi r)={T[0]*np.sqrt(2*np.pi*0.02):.4f}")
