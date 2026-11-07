import numpy as np, time, sys
from crack_solver2 import Params, solve
P = Params(0.2, 0.4, 0.89)
s = float(sys.argv[1])
for N in [int(v) for v in sys.argv[2:]]:
    t0 = time.time()
    o = solve(s, N, P)
    o2 = solve(s, N, P, kmax=2*(4.0*(2*N+1) + 40.0*np.sqrt(s) + 200.0))
    print(f"s={s:8.1e} N={N:4d} sK1={s*o['L1']['K']:.8f} (2kmax {s*o2['L1']['K']:.8f}) sV1={s*o['L1']['V']:.8f} sK2={s*o['L2']['K']:.8f} sV2={s*o['L2']['V']:.8f} |b_last|={abs(s*o['L1']['b'][-1]):.1e} ({time.time()-t0:.1f}s)", flush=True)
