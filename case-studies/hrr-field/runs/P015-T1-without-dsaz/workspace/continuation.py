"""
Continuation in n: start at n=1.5 (converged), step up to n=6.5 using the
previous solution as the seed. Only one axis BC: phi'(0)=0; free ICs
(phi''(0), phi'''(0)) solved from A(pi)=0, C(pi)=0.
"""
import numpy as np
import json
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import hrr2 as H


def endvals(n, u):
    H.n = n
    H.p = 1.0/(n+1.0)
    H.k = 2.0 - H.p
    Y0 = [1.0, 0.0, u[0], u[1]]
    s = solve_ivp(H.rhs, [0.0, np.pi], Y0, rtol=1e-10, atol=1e-12,
                  method='Radau', max_step=0.01)
    if not s.success:
        return np.array([1e3, 1e3])
    phi, d1, d2, d3 = s.y[:, -1]
    return np.array([H.k*phi + d2, -(H.k-1.0)*d1])


def solve(n, guess):
    sol, info, ier, msg = fsolve(lambda x: endvals(n, x), guess,
                                 xtol=1e-12, full_output=True)
    r = endvals(n, sol)
    return sol, np.linalg.norm(r), ier


ns = np.arange(1.5, 6.5001, 0.1)
guess = np.array([-1.125120, 0.653734])
ok = True
results = []
for n in ns:
    sol, res, ier = solve(n, guess)
    if res > 1e-6:
        ok = False
        print("FAIL at n=%.3f res=%.2e ier=%d" % (n, res, ier))
        break
    guess = sol
    results.append((n, sol[0], sol[1], res))
    if abs(n - round(n)) < 1e-9 and len(results) % 5 == 0:
        print("n=%.3f  phi''(0)=%.6f phi'''(0)=%.6f res=%.2e"
              % (n, sol[0], sol[1], res), flush=True)

# final report at n=6.5
sol, res, ier = solve(6.5, guess)
print("FINAL n=6.5  phi''(0)=%.10f phi'''(0)=%.10f res=%.2e ier=%d"
      % (sol[0], sol[1], res, ier))
np.save("phi_ic_6p5.npy", sol)
# also save the angular solution for later use
H.n = 6.5
H.p = 1.0/(6.5+1.0)
H.k = 2.0 - H.p
Y0 = [1.0, 0.0, sol[0], sol[1]]
s = solve_ivp(H.rhs, [0.0, np.pi], Y0, rtol=1e-12, atol=1e-14,
              method='Radau', dense_output=True, max_step=0.005)
np.savez("hrr_sol_6p5.npz", th=s.t, Y=s.y.T)
print("saved solution; points", s.t.size)