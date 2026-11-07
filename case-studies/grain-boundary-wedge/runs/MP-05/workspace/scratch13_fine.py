import numpy as np
from solver8 import Resolvent
import warnings; warnings.filterwarnings('ignore')

N = 28
R = Resolvent(N)
def smin(p):
    A = np.hstack([(p*R.F0 - 1.0)[:, None], p*R.F - R.H])
    return np.linalg.svd(A, compute_uv=False)[-1]

lam = np.logspace(np.log10(0.05), np.log10(5000), 25000)
vals = np.array([smin(l) for l in lam])
dips = []
for i in range(1, len(lam)-1):
    if vals[i] < vals[i-1] and vals[i] < vals[i+1] and vals[i] < 0.05*min(vals[i-1], vals[i+1]):
        dips.append(i)
print("dip count:", len(dips))
from scipy.optimize import minimize_scalar
found = []
for i in dips:
    lo, hi = lam[max(i-2,0)], lam[min(i+2, len(lam)-1)]
    r = minimize_scalar(lambda l: np.log10(smin(l)), bounds=(lo, hi), method='bounded',
                        options={'xatol': 1e-11})
    found.append((r.x, 10**r.fun))
for l, v in found:
    print(f"  lambda = {l:.8f}   smin = {v:.3e}")
