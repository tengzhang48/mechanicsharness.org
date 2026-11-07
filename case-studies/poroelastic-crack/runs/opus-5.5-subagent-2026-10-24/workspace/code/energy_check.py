"""Check of the Laplace-domain J-integral / energy-release identity for loading 1
(homogeneous Dirichlet pore-pressure condition on the faces):
   J_tip = (1-nu) pi Kn^2/2 + pi (sum b_n)^2/(2 kappa s)      (G=1, Kn = K/sqrt(pi L))
   J_energy = (1/(4s)) dV/dL = V/s + (1/2) dV/ds     (L = 1, scaling V = L^2 s^-1 F(s L^2))
"""
import numpy as np
from crack_solver2 import Params, solve
P = Params(0.2, 0.4, 0.89)
for s in [0.1, 1.0, 10.0, 100.0]:
    N = 30
    o = solve(s, N, P)['L1']
    h = 1e-4*s
    Vp = solve(s + h, N, P)['L1']['V']; Vm = solve(s - h, N, P)['L1']['V']
    dV = (Vp - Vm)/(2*h)
    Jen = o['V']/s + 0.5*dV
    Jmech = (1 - P.nu)*np.pi*o['K']**2/2
    Jfl = np.pi*np.sum(o['b'])**2/(2*P.kappa*s)
    print(f"s={s:7.2f}  s^2*J_energy={s*s*Jen:.8f}  s^2*(Jmech+Jfluid)={s*s*(Jmech+Jfl):.8f}  (mech {s*s*Jmech:.6f}, fluid {s*s*Jfl:.6f})")
