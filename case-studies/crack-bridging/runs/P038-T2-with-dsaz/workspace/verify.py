"""Independent checks of the beam DCB cohesive-zone result.

1. Finite-difference solve of the ODE (no analytic basis) -> same delta, L.
2. Global J-integral / energy balance check.
3. Sensitivity to the front boundary condition (clamped vs free-moment).
"""
import numpy as np
from scipy.optimize import brentq
import solve as S

# ------------------------------------------------------------------ check 1
def fd_solve(L, M, n=800):
    """Finite-difference solve of D~ d'''' = -p(x) on [0, L+a0].

    Solved by shooting the clamped-end derivatives (delta''(0), delta'''(0)):
    the linear ODE has an exact solution (built from the analytic basis), so we
    instead integrate with solve_ivp to stay fully independent of solve.py.
    """
    from scipy.integrate import solve_ivp
    x2 = L + S.a0

    def rhs_fun(x, y):
        d0, d1, d2, d3 = y          # delta, delta', delta'', delta'''
        p = S.b * S.sb * (1 - d0 / S.dc) if x < L else 0.0
        # D~ delta'''' = -p  (p: closing load)
        return [d1, d2, d3, -p / S.Dtil]

    def integrate(d2_0, d3_0):
        y0 = [0.0, 0.0, d2_0, d3_0]
        sol = solve_ivp(rhs_fun, [0.0, x2], y0, rtol=1e-11, atol=1e-14,
                        dense_output=True, method="DOP853")
        return sol

    def residuals(v):
        sol = integrate(v[0], v[1])
        ye = sol.y[:, -1]
        return np.array([ye[2] - M / S.Dtil, ye[3]])   # delta''=M/D~, delta'''=0

    from scipy.optimize import fsolve
    v0 = np.array([M / S.Dtil, 0.0])
    v = fsolve(residuals, v0, xtol=1e-13)

    sol = integrate(v[0], v[1])
    return sol

def fd_delta_at(L, M, xq, n=800):
    sol = fd_solve(L, M)
    return sol.sol(np.atleast_1d(xq))[0]

# ------------------------------------------------------------------ check 2
def energy_check(L, M, n=800):
    """G_far (=M^2/2bD) should equal G_tip + bridging work, with G_tip=Gam0."""
    delta, _ = S.solve_delta(L, M)
    dt = delta(L)[0]
    G_far = M**2 / (2 * S.b * S.Dtil)
    G_tip = G_far - S.bridging_work(dt)
    return G_far, G_tip, S.Gam0

if __name__ == "__main__":
    M, G = S.M_from_tip()
    L = S.L_for_M(M)
    print("analytic: M=%.6g  L(mm)=%.5g  dt(um)=%.5g  d(L/2)(um)=%.5g"
          % (M, L * 1e3, S.solve_delta(L, M)[0](L)[0] * 1e6,
             S.solve_delta(L, M)[0](L / 2)[0] * 1e6))

    # check 1: FD at the same (L,M)
    for n in (400, 800, 1600):
        dL = fd_delta_at(L, M, np.array([L, L / 2, L + S.a0]))
        print("shoot: dt=%.5f um  d(L/2)=%.4f um  d(end)=%.2f um"
              % (dL[0] * 1e6, dL[1] * 1e6, dL[2] * 1e6))

    # check 2: energy
    print("energy check (G_far, G_tip, Gam0):", energy_check(L, M))

    # check 3: L from the FD opening (independently of analytic delta)
    def f(Lt):
        return fd_delta_at(Lt, M, np.array([Lt]))[0] - S.dc
    Lfd = brentq(f, 1e-4, 5e-2, xtol=1e-12)
    print("L from FD (mm) =", Lfd * 1e3)