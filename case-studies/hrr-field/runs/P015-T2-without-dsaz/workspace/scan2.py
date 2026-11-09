import numpy as np
from scipy.optimize import brentq
import importlib.util
spec = importlib.util.spec_from_file_location("fv", "fast.py")
fv = importlib.util.module_from_spec(spec); spec.loader.exec_module(fv)

def endy(s, n, alpha, f2_0):
    sol = fv.integrate(s, n, alpha, f2_0)
    if (not sol.success) or sol.y.shape[1] == 0 or not np.all(np.isfinite(sol.y[:, -1])):
        return None
    return sol.y[:, -1]

def Rt(s, n, alpha, f2):
    y = endy(s, n, alpha, f2)
    return np.nan if y is None else (s-1.0)*y[1]

def Tt(s, n, alpha, f2):
    y = endy(s, n, alpha, f2)
    return np.nan if y is None else (2.0-s)*(1.0-s)*y[0]

def branches(s, n, alpha, lo=-10, hi=10, M=201):
    f2g = np.linspace(lo, hi, M)
    v = np.array([Rt(s, n, alpha, f) for f in f2g])
    out = []
    for i in range(M-1):
        if np.isfinite(v[i]) and np.isfinite(v[i+1]) and v[i]*v[i+1] < 0:
            try:
                r = brentq(lambda z: Rt(s, n, alpha, z), f2g[i], f2g[i+1], xtol=1e-13)
                out.append((r, Tt(s, n, alpha, r)))
            except Exception:
                pass
    return out

if __name__ == "__main__":
    n = 4.0
    alpha = 2*np.pi/3
    print(f"target alpha={alpha:.5f}, n={n}")
    for s in np.arange(0.05, 0.9, 0.025):
        b = branches(s, n, alpha)
        if b:
            txt = "; ".join(f"f2={r:.3f}:Tt={t:+.4f}" for r, t in b)
            print(f"s={s:.3f}: {txt}", flush=True)