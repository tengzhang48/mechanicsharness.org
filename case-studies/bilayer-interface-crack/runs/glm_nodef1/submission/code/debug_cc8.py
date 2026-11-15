import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import plane_strain_D, BilayerCrack
from central_crack import CentralCrack
from scipy.sparse.linalg import spsolve, eigsh

nu = 1/3
D = [plane_strain_D(1.0, nu)]*2

for name, mesh in [("CentralCrack", CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)),
                   ("BilayerCrack", BilayerCrack(1.0, 2.0, 30.0, 3e-4, 2e-4))]:
    K = fem.assemble(mesh.coords, mesh.elems, mesh.mats, D)
    A = mesh.nid[0, 0]
    pin = [2*A, 2*A+1, 2*mesh.nid[0, mesh.ny-1]]
    free = np.setdiff1d(np.arange(2*len(mesh.coords)), pin)
    Kf = K[free][:, free]
    f = np.zeros(2*len(mesh.coords))
    # probe: unit y-force at an interface node near the crack middle
    if name == "CentralCrack":
        pair = min(mesh.pairs, key=lambda p: abs(p[2]+1.0)); b = pair[0]
    else:
        b = mesh.nid[np.searchsorted(mesh.X, -1.0), mesh.j0]
    f[2*b+1] = 1.0
    u = np.zeros(2*len(mesh.coords))
    u[free] = spsolve(Kf.tocsc(), f[free])
    res = np.linalg.norm(K @ u - f)
    print(f"{name}: res={res:.2e}  max|u|={np.abs(u).max():.3e}  "
          f"nnz={Kf.nnz}  ndof={len(free)}")
