"""Scan the V1 stress-function eigenproblem g(pi; mu) = 0 over mu in (0.02, 2.5) for other symmetric
branches (audit 2026-11-07). mu < 1: singular stress r^(mu-1); finite nonzero J needs mu = n/(n+1)."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recompute_V1 as v1
mus = np.concatenate([np.linspace(0.02, 0.995, 196), np.linspace(1.005, 2.5, 300)])
vals = []
for m in mus:
    try:
        s = v1.solve_angular(m)[0]
        vals.append(float(s.y[0, -1]) if s.success and s.t[-1] >= np.pi - 1e-12 else None)
    except Exception:
        vals.append(None)
roots = [(float(mus[i]), float(mus[i + 1])) for i in range(len(mus) - 1)
         if vals[i] is not None and vals[i + 1] is not None and vals[i] * vals[i + 1] < 0]
out = {"mu_grid": [float(m) for m in mus], "g_pi": vals, "sign_change_brackets": roots,
       "failed_points": int(sum(v is None for v in vals)), "n": v1.n, "mu_key": v1.n / (v1.n + 1)}
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "scan_V1_branches.json"), "w"), indent=1)
print("sign-change brackets:", roots, " failed points:", out["failed_points"])
