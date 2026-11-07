import numpy as np, time
from crack_solver import Params, solve
P = Params(0.2, 0.4, 0.89)
for s in [1e2, 1e3]:
    for N in [20, 40, 80, 120]:
        for (width, npts) in [(0.5, 12), (0.25, 12)]:
            t0 = time.time()
            o = solve(s, N, P, width=width, npts=npts)
            print(f"s={s:7.0e} N={N:4d} w={width} sK1={s*o['L1']['K']:.7f} sV1={s*o['L1']['V']:.7f} sK2={s*o['L2']['K']:.7f} sV2={s*o['L2']['V']:.7f}  ({time.time()-t0:.1f}s)", flush=True)
