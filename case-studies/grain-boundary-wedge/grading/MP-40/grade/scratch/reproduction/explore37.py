exec(open('<round-runs>/runs/MP-40/explore36.py').read().split('def study')[0])

import numpy as np
import mpmath as mp
mp.mp.dps = 70

Q = 80
tg, wg = np.polynomial.legendre.leggauss(Q)
XG = [mp.mpf(0.5*(t+1.0)) for t in tg]
WQ = [mp.mpf(0.5*w) for w in wg]


def Nphi_mp(x, c):
    c3 = poly_deriv(c, 3)
    p3x = mp.mpf(peval(c3, x))
    tot1 = mp.mpf(0); tot2 = mp.mpf(0)
    for xi, wi in zip(XG, WQ):
        p3v = mp.mpf(peval(c3, xi))
        tot1 += wi*(p3v - p3x)/(x - xi)
        K2 = 1/(x+xi) + 2*xi*(x-xi)/(x+xi)**3
        tot2 += wi*p3v*K2
    return tot1 + p3x*mp.log(x/(1-x)) - tot2


def study(N):
    basis = make_basis(N)
    n = len(basis); m = n
    nodes = [mp.mpf(str(0.5*(1-np.cos(np.pi*(i+1)/(m+1))))) for i in range(m)]
    On = mp.matrix(m, n); Vn = mp.matrix(m, n)
    for i, xi in enumerate(nodes):
        for k, p in enumerate(basis):
            On[i, k] = Nphi_mp(xi, p)
            Vn[i, k] = mp.polyval(list(p[::-1]), xi)
    L = mp.inverse(Vn)*On
    ev = mp.eig(L, left=False, right=False)
    cand = sorted([(-float(complex(e).real), k) for k, e in enumerate(ev)
                   if abs(complex(e).imag) < 1e-12*max(1, abs(complex(e))) and -float(complex(e).real) > 1e-6])
    return basis, nodes, L, ev, cand


print("N   n   lam1        lam2        lam3", flush=True)
Ns = [8, 12, 16, 20, 24, 28, 32, 36]
res = {}
for N in Ns:
    basis, nodes, L, ev, cand = study(N)
    res[N] = (basis, nodes, L, ev, cand)
    v = [c[0] for c in cand]
    print("%2d %3d  %s" % (N, len(basis), "  ".join("%.6f" % w for w in v[:3])), flush=True)

lam1 = np.array([res[N][4][0][0] for N in Ns])
print("\nlam1:", list(zip(Ns, ["%.6f" % v for v in lam1])))
print("diffs:", ["%.5f" % v for v in np.diff(lam1)])