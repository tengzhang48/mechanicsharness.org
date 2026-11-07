"""Time-domain K(t), V(t) for both loadings by fixed-Talbot inversion of the Laplace-domain
finite-crack solution.  t in units L^2/c.  Results -> time_domain.json"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np, json, sys
from multiprocessing import Pool
from crack_solver2 import Params, solve
from talbot import talbot_nodes

P = Params(0.2, 0.4, 0.89)
M = int(sys.argv[1]) if len(sys.argv) > 1 else 12
times = [1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0, 1e3, 1e4]

def Nfor(s):
    return max(16, int(20 + 2.2*abs(s)**0.25))

def evalF(s):
    o = solve(s, Nfor(s), P)
    return [o['L1']['K'], o['L1']['V'], o['L2']['K'], o['L2']['V']]

if __name__ == '__main__':
    jobs = []
    for t in times:
        r, S, sig = talbot_nodes(t, M)
        jobs.append(complex(r)); jobs += list(S)
    with Pool(30) as pool:
        vals = pool.map(evalF, jobs)
    out = []; i = 0
    for t in times:
        r, S, sig = talbot_nodes(t, M)
        Fr = np.array(vals[i]); FS = np.array(vals[i+1:i+M]); i += M
        res = []
        for q in range(4):
            tot = 0.5*np.exp(r*t)*Fr[q] + np.sum(np.real(np.exp(t*S)*FS[:, q]*(1 + 1j*sig)))
            res.append(float(np.real((r/M)*tot)))
        out.append(dict(t=t, K1=res[0], V1=res[1], K2=res[2], V2=res[3]))
        print(json.dumps(out[-1]), flush=True)
    json.dump(out, open(f'time_domain_M{M}.json', 'w'), indent=1)
