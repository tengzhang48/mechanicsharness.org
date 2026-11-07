import numpy as np, time, sys, json
from crack_solver2 import Params, solve
P = Params(0.2, 0.4, 0.89)
s = float(sys.argv[1]); res = []
for N in [int(v) for v in sys.argv[2:]]:
    t0 = time.time()
    o = solve(s, N, P)
    r = dict(s=s, N=N, sK1=s*o['L1']['K'], sV1=s*o['L1']['V'], sK2=s*o['L2']['K'], sV2=s*o['L2']['V'],
             blast=float(abs(s*o['L1']['b'][-1])), time=time.time()-t0)
    res.append(r); print(json.dumps(r), flush=True)
json.dump(res, open(f'large_s_{s:.0e}.json', 'w'), indent=1)
