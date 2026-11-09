import numpy as np, time
from notch import integrate

def scan(n, alpha, s, pmin=-5.0, pmax=3.0, N=60):
    ps = np.linspace(pmin, pmax, N)
    vals = []
    for p in ps:
        try:
            vals.append(integrate(n, s, p, alpha).y[0, -1])
        except Exception:
            vals.append(np.nan)
    return ps, np.asarray(vals)

t0=time.time()
for s in [1.8]:
    ps, vals = scan(4.0, np.pi, s)
    print('crack dir:', [(round(ps[i],3), round(vals[i],3)) for i in range(len(ps))])
print('elapsed', time.time()-t0)
