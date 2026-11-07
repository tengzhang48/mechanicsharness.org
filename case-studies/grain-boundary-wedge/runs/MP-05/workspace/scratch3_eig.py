import numpy as np
from core import *

for N in [16, 24, 32, 48, 64]:
    lam, V, res, A, B, xc = solve_eig(N)
    ok = [(np.real(l), np.imag(l), r) for l, r in zip(lam, res)
          if abs(np.imag(l)) < 1e-9 * max(abs(np.real(l)), 1e-30) and np.real(l) > 1e-6 and r < 1e-8]
    ok.sort(key=lambda t: t[0])
    top = [(f"{l:.10f}", f"{r:.1e}") for l, _, r in ok[:5]]
    print(f"N={N:3d}: nphys={len(ok):3d} first: {top}")
