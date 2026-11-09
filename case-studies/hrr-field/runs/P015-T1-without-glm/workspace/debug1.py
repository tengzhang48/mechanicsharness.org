"""Debug: (a) validate F_phi4 at n=1 against elastic cos^3(theta/2);
(b) shoot at n=6.5: scan (phi''(0), p), find roots of [phi(pi), phi'(pi)]."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from hrrlib import F_phi4, PI, n_t

# ---- (a) n=1 elastic IVP test ----
n, p = 1.0, 0.5
def ode1(th, y):
    return np.array([y[1], y[2], y[3], F_phi4(y[0], y[1], y[2], y[3], n, p)])
sol1 = solve_ivp(ode1, [0, PI], [1.0, 0.0, -0.75, 0.0], rtol=1e-12, atol=1e-13,
                 dense_output=True, method='LSODA')
th = np.linspace(0, PI, 21)
exact = np.cos(th / 2) ** 3
print("n=1 test: max|phi - cos^3| =", np.max(np.abs(sol1.sol(th)[0] - exact)))
print("n=1 test: phi(pi) =", sol1.sol(PI)[0], " phi'(pi) =", sol1.sol(PI)[1])

# ---- (b) shooting at n=6.5 ----
def shoot(x, n=6.5, ret_full=False):
    s0_, pp = x
    def ode(th, y):
        return np.array([y[1], y[2], y[3], F_phi4(y[0], y[1], y[2], y[3], n, pp)])
    s = solve_ivp(ode, [0, PI], [1.0, 0.0, s0_, 0.0], rtol=1e-11, atol=1e-12,
                  dense_output=True, method='LSODA')
    yb = s.y[:, -1]
    if ret_full:
        return s
    return yb[0], yb[1]

# root find from scan best
from scipy.optimize import root
for start in [(-0.657, 0.16), (-0.657, 0.14), (-0.6, 0.1333)]:
    r = root(lambda x: shoot(x), start, method='hybr', tol=1e-13)
    print("root from", start, "->", r.success, r.x, shoot(r.x))

# coarse scan
best = []
for p0 in np.linspace(0.10, 0.17, 29):
    for s0_ in np.linspace(-1.6, 0.6, 29):
        try:
            ph, dph = shoot((s0_, p0))
        except Exception:
            continue
        best.append((abs(ph) + abs(dph), s0_, p0, ph, dph))
best.sort(key=lambda t: t[0])
print("best scan points (norm, phi''0, p, phi(pi), phi'(pi)):")
for b in best[:8]:
    print(f"  {b[0]:.4f}  phi''0={b[1]:+.3f} p={b[2]:.4f} phi(pi)={b[3]:+.4f} phi'(pi)={b[4]:+.4f}")
