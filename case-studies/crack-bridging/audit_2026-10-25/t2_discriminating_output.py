"""Which T2 output discriminates arm-stiffness errors? delta(L/2)/delta_c depends only on
gamma = Gamma0/(sigma_b delta_c); the rear rotation delta'(L) and delta at a fixed x scale with ell ~ C^(-1/4).
Own solver (same as recompute_P038.py T2_solve). Run with python -I."""
import math, json, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
b, E1, h1, E2, h2 = 25.0e-3, 110e9, 0.60e-3, 190e9, 0.90e-3
sb, dc, G0 = 0.20e6, 0.20e-3, 1.40
def run(D1, D2, G0=G0):
    C = 1/D1 + 1/D2; ell = (dc/(C*sb))**0.25; kap = math.sqrt(2*C*G0)*ell**2/dc
    f = lambda z, y: [y[1], y[2], y[3], y[0]-1.0]
    def res(p):
        s = solve_ivp(f, (0, p[0]), [0, 0, kap, p[1]], method="DOP853", rtol=1e-13, atol=1e-15)
        return [s.y[0, -1]-1, s.y[3, -1]]
    lam, A = fsolve(res, (1.7, 1.2), xtol=1e-14)
    s = solve_ivp(f, (0, lam), [0, 0, kap, A], method="DOP853", rtol=1e-13, atol=1e-15, dense_output=True)
    return {"L_mm": lam*ell*1e3, "rear_rotation_mrad": s.y[1, -1]*dc/ell*1e3,
            "delta_at_5mm_um": s.sol(5e-3/ell)[0]*dc*1e6, "delta_mid_um": s.sol(lam/2)[0]*dc*1e6}
D1, D2 = E1*h1**3/12, E2*h2**3/12
Dav = ((E1+E2)/2)*((h1+h2)/2)**3/12
out = {"key_model": run(D1, D2), "averaged_equal_arms": run(Dav, Dav), "only_upper_arm": run(D1, 1e30),
       "Gamma0_omitted": run(D1, D2, G0=0.0)}
json.dump(out, open("t2_discriminating_output.json", "w"), indent=1); print(json.dumps(out, indent=1))
