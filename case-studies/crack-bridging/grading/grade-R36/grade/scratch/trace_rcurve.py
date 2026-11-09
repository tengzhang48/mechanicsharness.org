"""Trace the R-curve by CONTINUATION in L (avoid spurious roots).
For each L, G_tip = Gamma_0 gives M; continue M from the previous L."""
import numpy as np
import solve as S

def trace(Lmax=40e-3, n=80):
    M_init = np.sqrt(2 * S.b * S.Dtil * S.Gam0)
    rows = []
    M_prev = M_init
    Ls = np.geomspace(1e-6, Lmax, n)
    for L in Ls:
        def g(M):
            d, _ = S.solve_delta(L, M)
            return M**2 / (2 * S.b * S.Dtil) - S.bridging_work(d(L)[0]) - S.Gam0
        # find root near M_prev by scanning outward
        Ms = np.linspace(max(M_prev * 0.8, 1e-8), M_prev * 3.0, 2000)
        vals = np.array([g(m) for m in Ms])
        roots = []
        for i in range(len(Ms) - 1):
            if vals[i] * vals[i + 1] < 0:
                from scipy.optimize import brentq
                roots.append(brentq(g, Ms[i], Ms[i + 1]))
        if not roots:
            continue
        # choose the root closest to M_prev (continuation)
        M = min(roots, key=lambda r: abs(r - M_prev))
        M_prev = M
        d, _ = S.solve_delta(L, M)
        dt = d(L)[0]
        rows.append((L, M, dt))
    return rows

if __name__ == "__main__":
    rows = trace()
    prev_dt = -1
    for L, M, dt in rows:
        flag = "  <-- dt crosses dc" if (prev_dt < S.dc <= dt) else ""
        print("L=%8.4f mm  M=%.6f Nm  dt=%10.4f um%s"
              % (L * 1e3, M, dt * 1e6, flag))
        prev_dt = dt