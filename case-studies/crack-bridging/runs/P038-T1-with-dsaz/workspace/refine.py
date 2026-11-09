import numpy as np
from solve import solve_finite, a_half, ell, sigma0, delta_c, Gamma0
an = a_half/ell
K0n = np.sqrt(Gamma0/(sigma0*delta_c))
E_=1.0; sig1=1.0
N=400
def state(Ln, N=N):
    c = an + Ln
    s0 = solve_finite(c, Ln, sigma0_=sig1, remote=0.0, N=N, E_=E_)
    s1 = solve_finite(c, Ln, sigma0_=0.0, remote=1.0, N=N, E_=E_)
    Sig = (K0n - s0['K_tip'])/s1['K_tip']
    sol = solve_finite(c, Ln, sigma0_=sig1, remote=Sig, N=N, E_=E_)
    th_a = np.arccos((c - Ln)/c)      # x = a (inner edge of zone)
    dt = sol['delta_at'](np.array([th_a]))[0]
    return Sig, dt
"""Definitive solve.

Condition at first connector failure (monotonic branch): the connector
adjacent to the original slit, at |x| = a (trailing edge of the bridging
zone), reaches delta = delta_c while the crack front at |x| = c = a + L is
at the intrinsic threshold K_tip = K0 = sqrt(E' Gamma0).  Uniform closure
sigma0 over a < |x| < c.

Normalized: lengths by ell = E' delta_c / sigma0, stress by sigma0, E'=1.
  alpha = a/ell,  K0n = sqrt(Gamma0/(sigma0 delta_c)).
For each Ln: find Sigma s.t. K_tip = K0n (linearity), read delta(root).
Bisect Ln on delta(root) = 1.
"""
import numpy as np
from solve import solve_finite, a_half, ell, sigma0, delta_c, Gamma0

an = a_half / ell
K0n = np.sqrt(Gamma0 / (sigma0 * delta_c))
E_ = 1.0


def state(Ln, N):
    c = an + Ln
    s0 = solve_finite(c, Ln, sigma0_=1.0, remote=0.0, N=N, E_=E_)
    s1 = solve_finite(c, Ln, sigma0_=0.0, remote=1.0, N=N, E_=E_)
    Sig = (K0n - s0['K_tip']) / s1['K_tip']
    sol = solve_finite(c, Ln, sigma0_=1.0, remote=Sig, N=N, E_=E_)
    th_a = np.arccos((c - Ln) / c)          # x = a (root / trailing edge)
    dt = sol['delta_at'](np.array([th_a]))[0]
    return Sig, dt


def bisect(N, iters=22):
    lo, hi = 0.40, 0.70
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        _, dt = state(mid, N)
        if dt < 1.0:
            lo = mid
        else:
            hi = mid
    Ln = 0.5 * (lo + hi)
    Sig, dt = state(Ln, N)
    return Ln, Sig, dt


for N in [400, 700, 1000]:
    Ln, Sig, dt = bisect(N)
    print(f"N={N:5d}: Ln*={Ln:.6f} L={Ln*ell*1e3:.5f} mm  "
          f"Sigma={Sig*160:.4f} MPa  delta(root)={dt:.6f}")
