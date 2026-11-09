import numpy as np
from scipy.optimize import brentq
import importlib.util
spec = importlib.util.spec_from_file_location("fv", "fast.py")
fv = importlib.util.module_from_spec(spec); spec.loader.exec_module(fv)

def endy(s, n, al, f2):
    sol = fv.integrate(s, n, al, f2)
    if (not sol.success) or sol.y.shape[1] == 0: return None
    y = sol.y[:, -1]
    return y if np.all(np.isfinite(y)) else None

def Rt(s, n, al, f2):
    y = endy(s, n, al, f2)
    return np.nan if y is None else (s-1)*y[1]
def Tt(s, n, al, f2):
    y = endy(s, n, al, f2)
    return np.nan if y is None else (2-s)*(1-s)*y[0]

def crossings(s, n, al, lo=-8, hi=8, M=161):
    g = np.linspace(lo, hi, M)
    v = np.array([Rt(s,n,al,f) for f in g])
    out=[]
    for i in range(M-1):
        if np.isfinite(v[i]) and np.isfinite(v[i+1]) and v[i]*v[i+1] < 0:
            try:
                r = brentq(lambda z: Rt(s,n,al,z), g[i], g[i+1], xtol=1e-12)
                out.append((r, Tt(s,n,al,r)))
            except Exception: pass
    return out

if __name__ == "__main__":
    n=4.0
    al=np.pi*0.999
    print("crack alpha=pi : looking for Tt=0 among Rt=0 branches")
    for s in np.arange(0.03,2.0,0.03):
        c = crossings(s,n,al)
        for r,t in c:
            if abs(t) < 0.5:
                print(f"  near: s={s:.3f} f2={r:.4f} Tt={t:+.4f}", flush=True)
    print("done", flush=True)
    # summarize sign of Tt on the "negative-f2 physical-ish" branch
    for s in np.arange(0.03,0.62,0.03):
        c = crossings(s,n,al)
        if c:
            print(f"  s={s:.3f}: " + "; ".join(f"f2={r:.3f}:Tt={t:+.3f}" for r,t in c), flush=True)