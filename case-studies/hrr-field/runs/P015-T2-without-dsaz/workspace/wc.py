import numpy as np
from scipy.optimize import brentq
import importlib.util
spec = importlib.util.spec_from_file_location('fv', 'fast.py')
fv = importlib.util.module_from_spec(spec); spec.loader.exec_module(fv)

n = 4.0; al = np.pi*0.999

def end(s, f2):
    sol = fv.integrate(s, n, al, f2)
    if (not sol.success) or sol.y.shape[1] == 0: return None
    y = sol.y[:, -1]
    return y if np.all(np.isfinite(y)) else None

def Rt(s, f2):
    y = end(s, f2); return np.nan if y is None else (s-1)*y[1]

def Tt(s, f2):
    y = end(s, f2); return np.nan if y is None else (2-s)*(1-s)*y[0]

g = np.linspace(-60, 60, 481)
for s in [0.1, 0.2, 0.25, 0.3]:
    v = np.array([Rt(s, f) for f in g])
    c = [(g[i], g[i+1]) for i in range(len(g)-1)
         if np.isfinite(v[i]) and np.isfinite(v[i+1]) and v[i]*v[i+1] < 0]
    print('s=', s, 'Rt-root branches:', len(c), flush=True)
    for lo, hi in c:
        try:
            r = brentq(lambda z: Rt(s, z), lo, hi)
            print('   f2=%.4f  Tt=%+.4f' % (r, Tt(s, r)), flush=True)
        except Exception as e:
            print('   err', e, flush=True)