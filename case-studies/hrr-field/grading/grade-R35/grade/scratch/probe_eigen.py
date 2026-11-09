import sys
import numpy as np
from solve import make_rhs, n, S_FORMULA
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def integ_fast(s, p, rtol=1e-9, max_step=0.01):
    rhs = make_rhs(n, s)
    f = lambda t, y: [y[1], y[2], y[3], rhs(y[0], y[1], y[2], y[3])]
    return solve_ivp(f, [0.0, np.pi], [1.0, 0.0, p, 0.0],
                     rtol=rtol, atol=1e-11, dense_output=False, max_step=max_step)

print("formula s =", S_FORMULA, flush=True)
print("--- f(pi) vs p at s=formula ---", flush=True)
for p in [-0.2, -0.4, -0.6, -0.8, -1.0, -1.2, -1.5, -2.0, -3.0, -5.0]:
    sol = integ_fast(S_FORMULA, p)
    print(f"p={p:6.2f}  f(pi)={sol.y[0,-1]: .6f}  fp(pi)={sol.y[1,-1]: .6f}", flush=True)

print("\n--- scan in s: find p with f(pi)=0, report f'(pi) ---", flush=True)
for ss in [1.3, 1.42, 1.5, 1.6, 1.7, 1.8, 1.8667, 1.95, 2.0]:
    ps = np.linspace(-6.0, -0.05, 40)
    vals = [integ_fast(ss, p).y[0, -1] for p in ps]
    roots = []
    for i in range(len(ps) - 1):
        if vals[i] * vals[i + 1] < 0:
            try:
                r = brentq(lambda p: integ_fast(ss, p).y[0, -1], ps[i], ps[i + 1],
                           xtol=1e-10)
                solr = integ_fast(ss, r)
                roots.append((r, solr.y[1, -1]))
            except Exception:
                pass
    print(f"s={ss:.4f}: " + ", ".join(f"(p={r:.4f}, fp(pi)={fp:.4f})"
                                       for r, fp in roots), flush=True)

print("\n--- 2D Newton near (1.8667, -0.697) ---", flush=True)
def bc(vec):
    ss, pp = vec
    sol = integ_fast(ss, pp, rtol=1e-11, max_step=0.005)
    return np.array([sol.y[0, -1], sol.y[1, -1]])

x = np.array([1.86666, -0.697])
for it in range(60):
    F = bc(x)
    J = np.zeros((2, 2))
    for j in range(2):
        dx = np.zeros(2); dx[j] = 1e-7
        J[:, j] = (bc(x + dx) - F) / 1e-7
    step = np.linalg.solve(J, -F)
    x = x + step
    if np.linalg.norm(step) < 1e-13:
        break
print(f"converged s={x[0]:.12f}, f''(0)={x[1]:.12f}, |F|={np.linalg.norm(F):.2e}, iters={it}", flush=True)
print(f"formula (2n+1)/(n+1) = {S_FORMULA:.12f}, diff = {abs(x[0]-S_FORMULA):.2e}", flush=True)