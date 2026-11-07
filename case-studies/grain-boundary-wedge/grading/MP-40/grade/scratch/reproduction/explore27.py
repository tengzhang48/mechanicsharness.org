exec(open('<round-runs>/runs/MP-40/explore26.py').read().split('prev = {}')[0])

import numpy as np
import mpmath as mp


def strong_residual(basis, coef, lam, xtest):
    # sigma=sum coef_k basis_k ; N sigma = sum coef_k N[basis_k]; check N sigma + lam sigma
    out = []
    for xv in xtest:
        xm = mp.mpf(xv)
        s = mp.mpf(0); ns = mp.mpf(0)
        for k, p in enumerate(basis):
            sk = mp.polyval(list(p[::-1]), xm)
            nk = mp.mpf(0)
            for d, c in enumerate(p):
                if abs(c) > 1e-14:
                    nk += mp.mpf(c)*Nval(d, xm)
            s += mp.mpf(float(coef[k]))*sk
            ns += mp.mpf(float(coef[k]))*nk
        out.append(float(ns+mp.mpf(lam)*s))
    return out


# ---- explore26 basis residual ----
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
cands = [(-complex(ev[k]).real, k) for k in range(n)
         if abs(complex(ev[k]).imag) < 1e-8*max(1, abs(complex(ev[k]))) and -complex(ev[k]).real > 1e-6]
cands.sort()
print("explore26 (w-basis) eigs:", ["%.5f" % c[0] for c in cands[:3]])
for lam, k in cands[:2]:
    coef = [W[k, j].real for j in range(n)]
    mx = max(abs(x) for x in coef)
    coef = [c/mx for c in coef]
    lamv = lam
    res = strong_residual(basis, coef, lamv, [0.35, 0.55, 0.75])
    scale = abs(mp.polyval(list(basis[0][::-1]), mp.mpf(0.5)))
    print("  lam=%.4f strong resid at .35/.55/.75: %s" % (lamv, ["%.3e" % r for r in res]))