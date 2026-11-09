"""Map R(phi''0, p) = (phi(pi), phi'(pi)); find sign crossings; 2D root."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from hrrlib import F_phi4, PI, n_t

def shoot(x, p, ret_full=False):
    def ode(th, y):
        return np.array([y[1], y[2], y[3], F_phi4(y[0], y[1], y[2], y[3], n_t, p)])
    s = solve_ivp(ode, [0, PI], [1.0, 0.0, x, 0.0], rtol=1e-11, atol=1e-13,
                  dense_output=True, method='LSODA')
    if ret_full:
        return s
    return s.y[0, -1], s.y[1, -1], s.status, np.max(np.abs(s.y[0]))

print("p      x1(phi_pi=0)   x2(phi'_pi=0)   x1-x2")
for p in np.linspace(0.08, 0.22, 29):
    xs = np.linspace(-2.0, 1.0, 121)
    vals = np.array([shoot(x, p)[:2] for x in xs])
    x1 = x2 = None
    for i in range(len(xs) - 1):
        if x1 is None and vals[i, 0] * vals[i + 1, 0] < 0:
            x1 = brentq(lambda x: shoot(x, p)[0], xs[i], xs[i + 1], xtol=1e-12)
        if x2 is None and vals[i, 1] * vals[i + 1, 1] < 0:
            x2 = brentq(lambda x: shoot(x, p)[1], xs[i], xs[i + 1], xtol=1e-12)
    if x1 is not None and x2 is not None:
        print(f"{p:.4f}  {x1:+.6f}   {x2:+.6f}   {x1 - x2:+.6f}")
    else:
        print(f"{p:.4f}  x1={x1}  x2={x2}")
