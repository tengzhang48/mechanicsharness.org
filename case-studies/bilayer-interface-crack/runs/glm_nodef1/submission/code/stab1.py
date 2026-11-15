"""Stability study: fit windows, mesh levels, dK extraction, solve residuals."""
import numpy as np
import sys
sys.path.insert(0, '/work/code')
from scipy.sparse.linalg import norm as spnorm
from fem import (BilayerCrack, plane_strain_D, end_loads, solve, assemble,
                 interface_tractions, crack_opening)
from fits import K_from_traction, amplitude_from_delta

nu = 1/3
E1, E2 = 9.0, 1.0
h, H = 1.0, 2.0
D = [plane_strain_D(E1, nu), plane_strain_D(E2, nu)]
Eb1, Eb2 = E1/(1-nu**2), E2/(1-nu**2)
S = 1/Eb1 + 1/Eb2
eps = np.log(2/3)/(2*np.pi)      # -0.0645033 (confirmed)
# delta relation: delta = Q * conj(K) * r^{1/2 - i eps}, Q from eigen-analysis
eps_m = -eps                     # +0.0645033
C = 4*S
Qpred = 1j*C/((1+2j*eps_m)*np.sqrt(2*np.pi)*np.cosh(np.pi*eps_m))
print(f"eps = {eps:.7f}, Qpred = {Qpred:.6f}")


def sigma(P, M):
    return {1: lambda x2: -P/h - 12*M*(x2 - h/2)/h**3,
            2: lambda x2:  P/H + 12*(M + P*(h+H)/2)*(x2 + H/2)/H**3}


def run(P, M, d1, d2, L=30.0, q1=1.35, q2=1.5):
    mesh = BilayerCrack(h, H, L, d1, d2, q1=q1, q2=q2, cap2=0.35)
    f = end_loads(mesh, sigma(P, M))
    u = solve(mesh, D, f)
    Kmat = assemble(mesh.coords, mesh.elems, mesh.mats, D)
    ru = Kmat @ u - f
    res = np.linalg.norm(ru)/max(np.linalg.norm(f), 1e-300)
    return mesh, u, res


print("case B (M=1):")
for (d1, d2) in [(1.2e-4, 3e-4), (6e-5, 1.5e-4), (3e-5, 7.5e-5),
                 (1.5e-5, 3.75e-5), (7.5e-6, 1.875e-5)]:
    mesh, u, res = run(0.0, 1.0, d1, d2)
    out = []
    for (rmin, rmax) in [(2e-3, 4e-2), (3e-3, 6e-2), (1e-3, 2e-2)]:
        rs = np.geomspace(rmin, rmax, 36)
        T = interface_tractions(mesh, D, u, rs)
        KT = K_from_traction(rs, T, eps, (0.0, 0.5, 1.0))
        rd, dd = crack_opening(mesh, u)
        m = (rd >= rmin) & (rd <= rmax)
        Ad = amplitude_from_delta(rd[m], dd[m], -eps, (0.0, 0.5, 1.0))
        Kd = np.conj(Ad/Qpred)
        out.append((rmin, rmax, KT, Kd))
    print(f" d1={d1:.1e} d2={d2:.1e} ndof={2*len(mesh.coords)} solver_res={res:.1e}")
    for (rmin, rmax, KT, Kd) in out:
        print(f"   win[{rmin:.0e},{rmax:.0e}]: KT={KT:.5f} ph={np.degrees(np.angle(KT)):+.3f} | "
              f"Kd={Kd:.5f} ph={np.degrees(np.angle(Kd)):+.3f}")
