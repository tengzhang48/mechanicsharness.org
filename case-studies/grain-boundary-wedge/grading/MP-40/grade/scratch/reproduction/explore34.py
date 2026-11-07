exec(open('<round-runs>/runs/MP-40/explore33.py').read().split('# cross-check')[0])

import numpy as np
import sympy as sp

X = sp.symbols('x')


def legendre_mono(j):
    e = sp.Poly(sp.expand(sp.legendre(j, 2*X-1)), X)
    return [float(e.coeff_monomial(X**d)) for d in range(e.degree()+1)]


xw3 = [0.0, 1.0, -3.0, 3.0, -1.0]


def make_basis(N):
    b = [[0.0, 3.0, -3.0, 1.0]]  # phi0 = 3x-3x^2+x^3
    for j in range(N+1):
        pj = legendre_mono(j)
        prod = np.convolve(xw3, pj)  # ascending
        b.append(list(prod))
    return b


def study(N):
    basis = make_basis(N)
    n = len(basis); m = n
    nodes = [0.5*(1-np.cos(np.pi*(i+1)/(m+1))) for i in range(m)]
    On = np.zeros((m, n)); Vn = np.zeros((m, n))
    for i in range(m):
        xi = nodes[i]
        for k, p in enumerate(basis):
            On[i, k] = Nphi_fast(xi, p)
            Vn[i, k] = peval(p, xi)
    L = np.linalg.solve(Vn, On)
    ev = np.linalg.eigvals(L)
    pos = sorted([-e.real for e in ev if abs(e.imag) < 1e-6*max(1, abs(e)) and -e.real > 1e-6])
    return basis, nodes, On, Vn, L, pos


print("N   n   lam1        lam2        lam3        lam4", flush=True)
Ns = [6, 9, 12, 15, 18, 21, 24, 27, 30, 34, 38]
res = {}
for N in Ns:
    basis, nodes, On, Vn, L, pos = study(N)
    res[N] = (basis, nodes, On, Vn, L, pos)
    print("%2d %3d  %s" % (N, len(basis), "  ".join("%.6f" % v for v in pos[:4])), flush=True)

lam1 = np.array([res[N][5][0] for N in Ns])
print("\nlam1 vs N:", list(zip(Ns, ["%.6f" % v for v in lam1])))
print("diffs:", ["%.5f" % v for v in np.diff(lam1)])