import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, solve, crack_opening
from central_crack import CentralCrack, end_loads_uniform
from scipy.sparse.linalg import spsolve
from scipy.sparse import identity

nu = 1/3
Db = [plane_strain_D(9.0, nu), plane_strain_D(1.0, nu)]

for (L, a, d1, d2) in [(10.0, 0.5, 1e-2, 2e-3), (10.0, 0.5, 1e-3, 5e-4),
                       (40.0, 1.0, 3e-4, 2e-4)]:
    mesh = CentralCrack(1.0, 2.0, L, d1, d2, a)
    f = end_loads_uniform(mesh, 1.0, both=True)
    K = fem.assemble(mesh.coords, mesh.elems, mesh.mats, Db)
    A = mesh.nid[0, 0]
    pin = [2*A, 2*A+1, 2*mesh.nid[0, mesh.ny-1]]
    free = np.setdiff1d(np.arange(2*len(mesh.coords)), pin)
    Kf = K[free][:, free].tocsc()
    for reg in [0.0, 1e-10]:
        Kr = (Kf + reg*identity(Kf.shape[0])).tocsc()
        u = np.zeros(2*len(mesh.coords))
        u[free] = spsolve(Kr, f[free])
        rd, dd = crack_opening(mesh, u)
        m = np.argmin(np.abs(rd-a*0.9))
        print(f"L={L} a={a} d1={d1:.0e} d2={d2:.0e} reg={reg:.0e}: "
              f"max|u|={np.abs(u).max():9.3f} delta(x=-{a*0.9:.2f})={dd[m]:+.3e}")
