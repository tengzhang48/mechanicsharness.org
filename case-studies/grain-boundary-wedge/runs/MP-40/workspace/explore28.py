exec(open('<round-runs>/runs/MP-40/explore26.py').read().split('prev = {}')[0])

import numpy as np
import mpmath as mp

M = 20
basis = make_basis(M)
n = len(basis); m = n
nodes = [mp.mpf(0.5)*(1-mp.cos(mp.pi*(i+1)/(m+1))) for i in range(m)]
On = mp.matrix(m, n); Vn = mp.matrix(m, n)
for i in range(m):
    xi = nodes[i]
    for k, p in enumerate(basis):
        acc = mp.mpf(0)
        for d, c in enumerate(p):
            if abs(c) > 1e-14:
                acc += mp.mpf(c)*Nval(d, xi)
        On[i, k] = acc
        Vn[i, k] = mp.polyval(list(p[::-1]), xi)
C = mp.inverse(Vn)*On
ev, W = mp.eig(C, left=False, right=True)
cands = sorted([(-complex(ev[k]).real, k) for k in range(n)
                if abs(complex(ev[k]).imag) < 1e-8*max(1, abs(complex(ev[k]))) and -complex(ev[k]).real > 1e-6])
lam, k1 = cands[0]
coef = [W[j, k1].real for j in range(n)]
mx = max(abs(c) for c in coef)
coef = [c/mx for c in coef]
def dot(row, i):
    return sum(row[i, j]*mp.mpf(coef[j]) for j in range(n))


print("lam=%.5f  residue at nodes:" % lam)
res = [float(dot(On, i) + mp.mpf(lam)*dot(Vn, i)) for i in range(m)]
print("  max node resid:", max(abs(r) for r in res))
# off-node ratio
for xv in [0.2, 0.35, 0.5, 0.65, 0.8, 0.9]:
    xm = mp.mpf(xv)
    s = mp.mpf(0); ns = mp.mpf(0)
    for k, p in enumerate(basis):
        nk = mp.mpf(0)
        for d, c in enumerate(p):
            if abs(c) > 1e-14:
                nk += mp.mpf(c)*Nval(d, xm)
        s += mp.mpf(coef[k])*mp.polyval(list(p[::-1]), xm)
        ns += mp.mpf(coef[k])*nk
    print("  x=%.2f  sigma=%.5f  Nsigma=%.5f  Nsigma/sigma=%.4f" % (xv, float(s), float(ns), float(ns/s) if abs(s) > 1e-9 else 0))

# also check node-values of sigma and Nsigma relationship to a TRUE eigenfunction
# print first few values
sig_nodes = [float(dot(Vn, i)) for i in range(m)]
print("sigma at first nodes:", ["%.4f" % v for v in sig_nodes[:6]])