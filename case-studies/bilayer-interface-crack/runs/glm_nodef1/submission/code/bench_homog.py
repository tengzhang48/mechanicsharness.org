"""Benchmarks of the FEM + extraction pipeline on homogeneous bilayers.

For a homogeneous strip with equal arms h = H and moments +-M (case B) the
exact result is pure mode I with K_I = 2 sqrt(3) M / h^{3/2}
(G = 12 M^2/(Ebar h^3) by the J-integral argument; K_I = sqrt(G Ebar)).
For case A (axial forces -P/+P) the exact G = P^2/(Ebar h) gives
|K| sqrt(h)/P = 1 with phase +-90 deg (pure mode II).
"""
import numpy as np
import sys
sys.path.insert(0, '/work/code')
from fem import (BilayerCrack, plane_strain_D, end_loads, solve,
                 interface_tractions, crack_opening)
from fits import measure_eps, K_from_traction, amplitude_from_delta

nu = 1/3
E = 1.0
h = H = 1.0
L = 30.0
D = [plane_strain_D(E, nu)]*2


def make_sigma(P, M, h, H):
    return {1: lambda x2: -P/h - 12*M*(x2 - h/2)/h**3,
            2: lambda x2:  P/H + 12*(M + P*(h+H)/2)*(x2 + H/2)/H**3}


def run(P, M, d1, d2):
    mesh = BilayerCrack(h, H, L, d1, d2)
    f = end_loads(mesh, make_sigma(P, M, h, H))
    u = solve(mesh, D, f)
    rs = np.geomspace(0.003, 0.12, 41)
    T = interface_tractions(mesh, D, u, rs)
    rd, dd = crack_opening(mesh, u)
    m = (rd > 0.003) & (rd < 0.12)
    eps, res, _ = measure_eps(rs, T)
    K = K_from_traction(rs, T, 0.0)
    A_d = amplitude_from_delta(rd[m], dd[m], 0.0)
    return eps, res, K, A_d, mesh, u


print("=== homogeneous, case B (M=1): exact K_I = 2*sqrt(3) = %.6f, phase 0" % (2*np.sqrt(3)))
for (d1, d2) in [(3e-4, 5e-3), (1.5e-4, 2.5e-3), (7.5e-5, 1.25e-3)]:
    eps, res, K, A_d, mesh, u = run(0.0, 1.0, d1, d2)
    print(f"  d1={d1:.1e} d2={d2:.1e} ndof={2*len(mesh.coords)}: "
          f"eps={eps:+.5f} |K|={abs(K):.5f} ph={np.degrees(np.angle(K)):+.3f} deg "
          f"res={res:.2e} A_d={A_d:.5f} (dof: |K|-exact)/exact = {abs(K)/(2*np.sqrt(3))-1:+.3%}")

print("=== homogeneous, case A (P=1): exact |K| = 1, phase +-90")
for (d1, d2) in [(3e-4, 5e-3), (1.5e-4, 2.5e-3)]:
    eps, res, K, A_d, mesh, u = run(1.0, 0.0, d1, d2)
    print(f"  d1={d1:.1e}: eps={eps:+.5f} |K|={abs(K):.5f} ph={np.degrees(np.angle(K)):+.3f} deg res={res:.2e}")
