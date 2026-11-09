import numpy as np
from solve import solve_finite, a_half, ell, sigma0, delta_c, Gamma0
an = a_half/ell
K0n = np.sqrt(Gamma0/(sigma0*delta_c))
E_=1.0; sig1=1.0
for Ln in [0.4, 0.5, 0.6]:
    c = an + Ln
    s0 = solve_finite(c, Ln, sigma0_=sig1, remote=0.0, N=400, E_=E_)
    s1 = solve_finite(c, Ln, sigma0_=0.0, remote=1.0, N=400, E_=E_)
    Sig = (K0n - s0['K_tip'])/s1['K_tip']
    sol = solve_finite(c, Ln, sigma0_=sig1, remote=Sig, N=400, E_=E_)
    xs = np.linspace(an, c, 11)
    th = np.arccos(xs/c)
    d = sol['delta_at'](th)
    print(f"Ln={Ln} Sig={Sig:.4f}")
    for xi,di in zip(xs,d):
        tag='ZONE' if xi>=an-1e-9 else 'crack'
        print(f"   x={xi:.4f}  delta={di:.4f}")
