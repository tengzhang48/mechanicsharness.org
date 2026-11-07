exec(open('<round-runs>/runs/MP-40/explore37.py').read().split('\ndef study')[0])

import numpy as np
import mpmath as mp


def study2(N):
    basis = make_basis(N)
    n = len(basis); m = n
    nodes = [mp.mpf(str(0.5*(1-np.cos(np.pi*(i+1)/(m+1))))) for i in range(m)]
    On = mp.matrix(m, n); Vn = mp.matrix(m, n)
    for i, xi in enumerate(nodes):
        for k, p in enumerate(basis):
            On[i, k] = Nphi_mp(xi, p)
            Vn[i, k] = mp.polyval(list(p[::-1]), xi)
    # column scaling (similarity: preserves eigenvalues)
    d = []
    for k in range(n):
        nk = mp.sqrt(sum(On[i, k]**2 for i in range(m)))
        d.append(mp.mpf(1)/nk if nk > 0 else mp.mpf(1))
        for i in range(m):
            On[i, k] *= d[k]; Vn[i, k] *= d[k]
    L = mp.inverse(Vn)*On
    ev = mp.eig(L, left=False, right=False)
    cand = sorted([(-float(complex(e).real), k) for k, e in enumerate(ev)
                   if abs(complex(e).imag) < 1e-12*max(1, abs(complex(e))) and -float(complex(e).real) > 1e-6])
    return basis, nodes, L, ev, cand, d


print("N   n   lam1        lam2        lam3", flush=True)
Ns = [8, 12, 16, 20, 24, 28, 32, 36, 40]
res = {}
for N in Ns:
    basis, nodes, L, ev, cand, d = study2(N)
    res[N] = (basis, nodes, L, ev, cand, d)
    print("%2d %3d  %s" % (N, len(basis), "  ".join("%.6f" % c[0] for c in cand[:8])), flush=True)

lam1 = np.array([res[N][4][0][0] for N in Ns])
print("\nlam1:", list(zip(Ns, ["%.6f" % v for v in lam1])))
print("diffs:", ["%.5f" % v for v in np.diff(lam1)])