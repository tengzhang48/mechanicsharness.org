import numpy as np
from gbcore import *

# V1: numeric A[1](x) vs closed form, incl. x -> 0+ limit -4
for N in [40, 80]:
    xg = cheb_lobatto_01(N)
    xc = xg[1:]
    Q = build_Q(xc, xg)
    num = Q @ np.ones(N + 1)
    ex = A_exact_const(xc)
    rel = np.abs(num - ex) / np.maximum(np.abs(ex), 1e-3)
    i = np.argmax(rel)
    print(f"N={N}: max rel err A[1] = {rel.max():.3e} at x={xc[i]:.4f}")
    print(f"   near-0 row x={xc[0]:.2e}: num={num[0]:.8f} exact={ex[0]:.8f} (limit -4)")
