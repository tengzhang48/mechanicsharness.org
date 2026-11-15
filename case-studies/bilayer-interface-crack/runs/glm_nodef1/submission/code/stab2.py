"""Careful convergence study for |K| and phase, cases A and B."""
import numpy as np
import sys, json
sys.path.insert(0, '/work/code')
from fem import (BilayerCrack, plane_strain_D, end_loads, solve,
                 interface_tractions, crack_opening)
from fits import fit_series

nu = 1/3
E1, E2 = 9.0, 1.0
h, H = 1.0, 2.0
D = [plane_strain_D(E1, nu), plane_strain_D(E2, nu)]
Eb1, Eb2 = E1/(1-nu**2), E2/(1-nu**2)
S = 1/Eb1 + 1/Eb2
eps = np.log(2/3)/(2*np.pi)
G = {"A": 1/(2*Eb1*h) + 1/(2*Eb2*H) + 6*((h+H)/2)**2/(Eb2*H**3),
     "B": 6*(1/(Eb1*h**3) + 1/(Eb2*H**3))}
print(f"eps={eps:.7f} G_A={G['A']:.7f} G_B={G['B']:.7f}")


def sigma(P, M):
    return {1: lambda x2: -P/h - 12*M*(x2 - h/2)/h**3,
            2: lambda x2:  P/H + 12*(M + P*(h+H)/2)*(x2 + H/2)/H**3}


def wfit(r, y, powers, w):
    A = np.stack([(r**p)*w for p in powers], axis=1)
    c, *_ = np.linalg.lstsq(A, y*w, rcond=None)
    return c[0]


def run(P, M, d1, d2):
    mesh = BilayerCrack(h, H, 30.0, d1, d2, q1=1.35, q2=1.5, cap2=0.35)
    f = end_loads(mesh, sigma(P, M))
    u = solve(mesh, D, f)
    rs = np.geomspace(4e-4, 8e-2, 60)
    T = interface_tractions(mesh, D, u, rs)
    rd, dd = crack_opening(mesh, u)
    m = (rd >= 4e-4) & (rd <= 8e-2)
    return rs, T, rd[m], dd[m], 2*len(mesh.coords)


levels = [(1.2e-4, 3.0e-4), (6e-5, 1.5e-4), (3e-5, 7.5e-5),
          (1.5e-5, 3.75e-5), (7.5e-6, 1.875e-5), (4e-6, 1e-5)]

store = {}
for case, (P, M) in [("A", (1.0, 0.0)), ("B", (0.0, 1.0))]:
    print(f"\n===== case {case}")
    seq = {}
    for (d1, d2) in levels:
        rs, T, rd, dd, ndof = run(P, M, d1, d2)
        y = T*np.sqrt(2*np.pi*rs)*np.exp(-1j*eps*np.log(rs))
        yd = dd/np.sqrt(rd)*np.exp(+1j*eps*np.log(rd))   # delta exp = -eps
        row = {}
        for (wname, w) in [("unw", np.ones_like(rs)), ("w", 1/np.sqrt(rs))]:
            for (rmin, rmax) in [(5e-4, 1e-2), (1e-3, 2e-2), (2e-3, 4e-2), (3e-3, 6e-2)]:
                m = (rs >= rmin) & (rs <= rmax)
                for nt, pw in [(3, (0.0, 0.5, 1.0)), (4, (0.0, 0.5, 1.0, 1.5))]:
                    K = wfit(rs[m], y[m], pw, w[m])
                    row[f"T{wname}{rmin:.0e}-{rmax:.0e}.{nt}"] = K
                md = (rd >= rmin) & (rd <= rmax)
                if md.sum() > 8:
                    wd = 1/np.sqrt(rd[md])
                    Ad = wfit(rd[md], yd[md], (0.0, 0.5, 1.0), wd)
                    Qpred = 1j*4*S/((1-2j*eps)*np.sqrt(2*np.pi)*np.cosh(np.pi*(-eps)))
                    row[f"D{rmin:.0e}-{rmax:.0e}"] = np.conj(Ad/Qpred)
        seq[(d1, d2)] = row
        print(f" d={d1:.1e}/{d2:.1e} ndof={ndof}: "
              f"Tunw5e-4-1e-2.3={row['Tunw5e-04-1e-02.3']:.6f} "
              f"Tw1e-3-2e-2.3={row['Tw1e-03-2e-02.3']:.6f} "
              f"Tunw2e-3-4e-2.3={row['Tunw2e-03-4e-02.3']:.6f} "
              f"D1e-3-2e-2={row.get('D1e-03-2e-02', float('nan')):.6f}")
    store[case] = {f"{a:.1e}/{b:.1e}": v for (a, b), v in seq.items()}

# Richardson on a representative estimator per case
print("\nRichardson extrapolation (Tunw 5e-4-1e-2, 3 terms):")
for case in ["A", "B"]:
    vals = [store[case][f"{a:.1e}/{b:.1e}"]["Tunw5e-04-1e-02.3"] for (a, b) in levels]
    for v in vals:
        print(f"   {case}: {v:.6f}  |K|={abs(v):.6f} ph={np.degrees(np.angle(v)):+.4f}")
    e2 = vals[-3]-vals[-4]; e1 = vals[-2]-vals[-3]; e0 = vals[-1]-vals[-2]
    print(f"   diffs: {e2:.2e} {e1:.2e} {e0:.2e}")

json.dump({c: {k: {kk: complex(vv) for kk, vv in r.items()} for k, r in v.items()}
           for c, v in store.items()}, open('/work/code/stab2_out.json', 'w'))
print("saved stab2_out.json")
