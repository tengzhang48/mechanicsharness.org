"""Long-time (small s) behaviour: s*K and s*V for both loadings at s -> 0,
fitted to  X(s) = X_inf + C/(ln(1/s) + D)  (2D diffusion: logarithmic approach)."""
import numpy as np
from crack_solver2 import Params, solve
P = Params(0.2, 0.4, 0.89)
ss = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12]
rows = []
for s in ss:
    o = solve(s, 16, P)
    rows.append([s, s*o['L1']['K'], s*o['L1']['V'], s*o['L2']['K'], s*o['L2']['V']])
    print("s=%.0e  sK1=%.9f  sV1=%.9f  sK2=%.9f  sV2=%.9f" % tuple(rows[-1]))
rows = np.array(rows)
eta = P.eta
print('expected: K1 -> 1, V1 -> pi(1-nu) = %.9f, K2 -> -eta = %.9f, V2 -> -pi*eta*(1-nu) = %.9f' % (np.pi*(1-P.nu), -eta, -np.pi*eta*(1-P.nu)))
for col, name in [(3, 'K2'), (4, 'V2')]:
    # fit X_inf, C, D using last three points
    from scipy.optimize import fsolve
    x = np.log(1/rows[-3:, 0]); y = rows[-3:, col]
    f = lambda v: v[0] + v[1]/(x + v[2]) - y
    sol = fsolve(f, [y[-1], 1.0, 0.0])
    print(name, 'fit X_inf=%.6f C=%.4f D=%.4f' % tuple(sol))
