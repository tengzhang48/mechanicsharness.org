"""Tip factor a(gamma) = lim_{s->inf} W(s), W = s*(K1 + 2 beta K2), and the tip volume coefficient
cW(gamma) = lim (s V_W - V_u) sqrt(s)/(1-nu_u), from the finite-crack solver; nu fixed at 0.2,
nu_u chosen to give gamma; alpha irrelevant (set 0.9)."""
import numpy as np, json, sys
from crack_solver2 import Params, solve
nu = 0.2
res = []
for gam in [float(g) for g in sys.argv[1:]]:
    nuu = nu + gam*(1 - nu)
    P = Params(nu, nuu, 0.9)
    row = dict(gamma=gam, data=[])
    for s, N in [(1e3, 40), (1e4, 40), (1e5, 60)]:
        o = solve(s, N, P)
        W = s*o['L1']['K'] + 2*P.beta*s*o['L2']['K']
        VW = s*o['L1']['V'] + 2*P.beta*s*o['L2']['V']
        row['data'].append(dict(s=s, W=W, cW=(VW - np.pi*(1 - nuu))*np.sqrt(s)/(1 - nuu),
                                sK2=s*o['L2']['K'], sV2=s*o['L2']['V']))
    W = [d['W'] for d in row['data']]
    # Richardson with s^-1/2 steps of factor sqrt(10)
    r = 10**-0.5
    row['a_rich'] = W[2] + (W[2] - W[1])*r/(1 - r)
    res.append(row); print(json.dumps(row), flush=True)
json.dump(res, open('a_of_gamma_%s.json' % '_'.join(sys.argv[1:]), 'w'), indent=1)
