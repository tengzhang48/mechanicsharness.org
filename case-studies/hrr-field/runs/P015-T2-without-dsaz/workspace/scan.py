import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import importlib.util
spec = importlib.util.spec_from_file_location("fv", "fast.py")
fv = importlib.util.module_from_spec(spec)
# fast.py has a __main__ guard; safe to exec
spec.loader.exec_module(fv)

# For fixed s, find f2_0 such that Sigma_rt(alpha)=0 ; then scan Tt(alpha) vs s.
def endstate(s, n, alpha, f2_0):
    sol = fv.integrate(s, n, alpha, f2_0)
    if not sol.success or sol.y.shape[1] == 0:
        return None
    return sol.y[:, -1]

def Rt(s, n, alpha, f2_0):
    y = endstate(s, n, alpha, f2_0)
    if y is None: return np.nan
    f0v, f1v, f2v, f3v, Vv = y
    return (s - 1.0)*f1v

def Tt(s, n, alpha, f2_0):
    y = endstate(s, n, alpha, f2_0)
    if y is None: return np.nan
    f0v, f1v, f2v, f3v, Vv = y
    return (2.0 - s)*(1.0 - s)*f0v

def scan(s, n, alpha, f2lo=-8, f2hi=8, M=81):
    f2g = np.linspace(f2lo, f2hi, M)
    vals = np.array([Rt(s, n, alpha, f) for f in f2g])
    # find sign changes -> roots of Rt in f2_0
    roots = []
    for i in range(M-1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i+1]) and vals[i]*vals[i+1] < 0:
            try:
                r = brentq(lambda z: Rt(s, n, alpha, z), f2g[i], f2g[i+1], xtol=1e-12)
                roots.append(r)
            except Exception:
                pass
    out = []
    for r in roots:
        out.append((r, Tt(s, n, alpha, r)))
    return out

if __name__ == "__main__":
    n = 4.0
    for alpha, tag in [(np.pi*0.999, "crack (want ~0.2)"), (2*np.pi/3, "target")]:
        print(f"=== alpha={alpha:.4f} {tag} ===")
        for s in np.linspace(0.05, 0.9, 18):
            res = scan(s, n, alpha)
            if res:
                print(f"  s={s:.3f}: " + ", ".join(f"f2*={r:.3f}->Tt={t:.4f}" for r,t in res))