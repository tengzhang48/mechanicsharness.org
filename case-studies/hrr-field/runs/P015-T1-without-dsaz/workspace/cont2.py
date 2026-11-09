"""
Faster continuation with per-step printing, looser IVP tolerances.
Only 1 axis BC (phi'(0)=0); free ICs (phi''(0), phi'''(0)).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import hrr2 as H


def endvals(n, u):
    H.n = n
    H.p = 1.0/(n+1.0)
    H.k = 2.0 - H.p
    Y0 = [1.0, 0.0, u[0], u[1]]
    s = solve_ivp(H.rhs, [0.0, np.pi], Y0, rtol=1e-8, atol=1e-10,
                  method='Radau', max_step=0.03)
    if not s.success:
        return np.array([1e6, 1e6])
    phi, d1, d2, d3 = s.y[:, -1]
    return np.array([H.k*phi + d2, -(H.k-1.0)*d1])


guess = np.array([-1.125120, 0.653734])   # n=1.5 solution
n = 1.5
out = []
while n <= 6.5 + 1e-9:
    sol, info, ier, msg = fsolve(lambda x: endvals(n, x), guess, xtol=1e-10,
                                 full_output=True)
    res = np.linalg.norm(endvals(n, sol))
    print("n=%.3f  phi''(0)=%.6f phi'''(0)=%.6f res=%.2e ier=%d"
          % (n, sol[0], sol[1], res, ier), flush=True)
    if res > 1e-5:
        print("  !! residual large; stopping", flush=True)
        break
    guess = sol
    out.append((n, sol[0], sol[1]))
    n = round(n + 0.1, 3)

np.save("cont2_out.npy", np.array(out))
if abs(out[-1][0] - 6.5) < 1e-6:
    H.n = 6.5; H.p = 1.0/(7.5); H.k = 2.0-H.p
    Y0 = [1.0, 0.0, guess[0], guess[1]]
    s = solve_ivp(H.rhs, [0.0, np.pi], Y0, rtol=1e-12, atol=1e-14,
                  method='Radau', dense_output=True, max_step=0.005)
    np.savez("hrr_sol_6p5.npz", th=s.t, Y=s.y.T)
    print("SAVED hrr_sol_6p5.npz points", s.t.size, flush=True)