exec(open('<round-runs>/runs/MP-40/explore26.py').read().split('prev = {}')[0])

import numpy as np
import mpmath as mp


def build(M):
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
    return basis, nodes, On, Vn


for M in [16, 20, 24]:
    basis, nodes, On, Vn = build(M)
    n = len(basis)
    L = mp.inverse(Vn)*On
    ev, ER = mp.eig(L)
    cands = sorted([(complex(ev[k]).real, k) for k in range(n)
                    if abs(complex(ev[k]).imag) < 1e-9*max(1, abs(complex(ev[k])))])
    # slowest decay = eigenvalue closest to zero (smallest |Re|)
    lam, k1 = min(cands, key=lambda c: abs(c[0]))
    vec = [ER[j, k1] for j in range(n)]
    mx = max(abs(v) for v in vec)
    vec = [v/mx for v in vec]

    def sig(xm):
        return sum(vec[k]*mp.polyval(list(basis[k][::-1]), xm) for k in range(n))

    def Nsig(xm):
        tot = mp.mpf(0)
        for k, p in enumerate(basis):
            acc = mp.mpf(0)
            for d, c in enumerate(p):
                if abs(c) > 1e-14:
                    acc += mp.mpf(c)*Nval(d, xm)
            tot += vec[k]*acc
        return tot

    print("M=%d  eig=%.5f   (ratio should be constant = eig)" % (M, lam), flush=True)
    for xv in [0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.97]:
        xm = mp.mpf(str(xv))
        s = sig(xm); ns = Nsig(xm)
        ratio = complex(ns/s) if abs(s) > 1e-30 else 0
        print("   x=%.2f  sigma=%.5f  Nsigma=%.5f  ratio=%.4f%+.4fj" % (xv, float(complex(s).real), float(complex(ns).real), ratio.real, ratio.imag), flush=True)