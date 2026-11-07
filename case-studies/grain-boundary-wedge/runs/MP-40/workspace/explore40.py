exec(open('<round-runs>/runs/MP-40/explore37.py').read().split('\ndef study')[0])

import numpy as np
import mpmath as mp

N = 20
basis = make_basis(N)
n = len(basis); m = n
nodes = [mp.mpf(str(0.5*(1-np.cos(np.pi*(i+1)/(m+1))))) for i in range(m)]
On = mp.matrix(m, n); Vn = mp.matrix(m, n)
for i, xi in enumerate(nodes):
    for k, p in enumerate(basis):
        On[i, k] = Nphi_mp(xi, p)
        Vn[i, k] = mp.polyval(list(p[::-1]), xi)
L = mp.inverse(Vn)*On
ev, ER = mp.eig(L)
cand = sorted([(-float(complex(e).real), k) for k, e in enumerate(ev)
               if abs(complex(e).imag) < 1e-12*max(1, abs(complex(e))) and -float(complex(e).real) > 1e-6])
lam1, k1 = cand[0]
print("lambda1 = %.6f" % lam1)
vec = [ER[j, k1] for j in range(n)]
mx = max(abs(v) for v in vec)
vec = [v/mx for v in vec]


def sig(xv):
    return sum(vec[k]*mp.polyval(list(basis[k][::-1]), mp.mpf(str(xv))) for k in range(n))


def re(z):
    return float(complex(z).real)


f1 = sig(1.0)
print("f(1)=%.6f" % re(f1))
for xv in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
    print("  f(%.2f)=%.6f   norm=%.6f" % (xv, re(sig(xv)), re(sig(xv)/f1)))
print("mode1_shape_half = %.6f" % re(sig(0.5)/f1))