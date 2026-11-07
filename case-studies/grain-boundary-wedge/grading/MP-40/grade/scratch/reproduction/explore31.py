exec(open('<round-runs>/runs/MP-40/explore26.py').read().split('prev = {}')[0])

import numpy as np
import mpmath as mp
import sympy as sp

X = sp.symbols('x')


def legendre_mono(j):
    # monomial coefficients (ascending) of P_j(2x-1)
    e = sp.Poly(sp.expand(sp.legendre(j, 2*X - 1)), X)
    return [float(e.coeff_monomial(X**d)) for d in range(e.degree() + 1)]


def make_basis(N):
    phi0 = [0, 3, -3, 1]  # x(1+(1-x)+(1-x)^2) = 3x-3x^2+x^3
    b = [phi0]
    xw3 = [0, 1, -3, 3, -1]  # x(1-x)^3
    for j in range(N + 1):
        pj = legendre_mono(j)
        prod = np.polymul(xw3, pj)  # descending? np.polymul uses descending
        prod = list(prod[::-1])  # to ascending
        b.append(prod)
    return b


def run(N, nodes_shift=1.0):
    basis = make_basis(N)
    n = len(basis); m = n
    nodes = [mp.mpf(0.5)*(1-mp.cos(mp.pi*(i+1)/(m+1))) for i in range(m)]
    # slight perturbation to avoid any accidental node-poly cancellation
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
    L = mp.inverse(Vn)*On
    ev, ER = mp.eig(L)
    cands = sorted([(-complex(ev[k]).real, k) for k in range(n)
                    if abs(complex(ev[k]).imag) < 1e-8*max(1, abs(complex(ev[k]))) and -complex(ev[k]).real > 1e-6])
    lams = [c[0] for c in cands]
    return basis, nodes, On, Vn, L, ev, ER, cands, lams


print("N   n   real-decay eigs (first 5)", flush=True)
store = {}
for N in [6, 10, 14, 18, 22, 26, 30]:
    basis, nodes, On, Vn, L, ev, ER, cands, lams = run(N)
    store[N] = (basis, nodes, On, Vn, ER, cands)
    print("%2d %3d  %s" % (N, len(basis), " ".join("%.5f" % v for v in lams[:5])), flush=True)

# mode shape from largest N
N = 30
basis, nodes, On, Vn, ER, cands = store[N]
n = len(basis)
lam1, k1 = cands[0]
vec = [ER[j, k1] for j in range(n)]
mx = max(abs(v) for v in vec)
vec = [v/mx for v in vec]


def sig(xv):
    return sum(vec[k]*mp.polyval(list(basis[k][::-1]), xv) for k in range(n))


f1 = sig(mp.mpf(1))
print("lam1=%.6f  f(1)=%.5f" % (lam1, float(f1)))
for xv in [0.0, 0.25, 0.5, 0.75, 0.9, 1.0]:
    print("  f(%.2f)=%.6f" % (xv, float(sig(mp.mpf(xv))/f1)))

# extrapolation lambda1 ~ L + c/N^p
Ns = sorted(store.keys())
Lv = np.array([abs(store[N][5][0][0]) for N in Ns])

# Aitken on last three
def aitken(a, b, c):
    return a - (b - a)**2/(c - 2*b + a)
print("raw lam1:", list(zip(Ns, ["%.6f" % v for v in Lv])))
# Richardson using p from first differences
d = np.diff(Lv)
print("diffs:", ["%.5f" % v for v in d])
print("aitken(last3):", aitken(*Lv[-3:]))