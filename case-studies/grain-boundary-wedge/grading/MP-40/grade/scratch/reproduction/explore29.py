exec(open('<round-runs>/runs/MP-40/explore26.py').read().split('prev = {}')[0])

import numpy as np
import mpmath as mp

M = 24
basis = make_basis(M)
n = len(basis); m = n
nodes = np.array([float(mp.mpf(0.5)*(1-mp.cos(mp.pi*(i+1)/(m+1)))) for i in range(m)])
On = np.zeros((m, n)); Vn = np.zeros((m, n))
for i in range(m):
    xi = nodes[i]
    for k, p in enumerate(basis):
        acc = 0.0
        for d, c in enumerate(p):
            if abs(c) > 1e-14:
                acc += c*float(Nval(d, mp.mpf(xi)))
        On[i, k] = acc
        Vn[i, k] = np.polyval(p[::-1], xi)
ev, W = np.linalg.eig(np.linalg.solve(Vn, On))
cands = sorted([(-ev[k].real, k) for k in range(n)
                if abs(ev[k].imag) < 1e-6*max(1, abs(ev[k])) and -ev[k].real > 1e-6])
print("eigs:", ["%.5f" % c[0] for c in cands[:4]])
lam, k1 = cands[0]
v = np.real(W[:, k1]); v = v/np.max(np.abs(v))


def sig(x):
    return sum(v[k]*np.polyval(basis[k][::-1], x) for k in range(n))


def Nsig(x):
    tot = 0.0
    for k, p in enumerate(basis):
        acc = 0.0
        for d, c in enumerate(p):
            if abs(c) > 1e-14:
                acc += c*float(Nval(d, mp.mpf(x)))
        tot += v[k]*acc
    return tot


print("x     sigma     Nsigma    ratio")
for x in [0.15, 0.3, 0.45, 0.6, 0.75, 0.9]:
    s = sig(x); ns = Nsig(x)
    print("  %.2f  %.5f  %.5f  %.4f" % (x, s, float(ns), float(ns)/s if abs(s) > 1e-12 else 0))
print("f(0.5)/f(1)=%.6f" % (sig(0.5)/sig(1.0)))
print("node residual max:", np.max(np.abs(On @ v + lam*(Vn @ v))))