import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss

x = sp.symbols('x', positive=True)
Q = 200
tg, wg = leggauss(Q)
xg = 0.5*(tg+1.0); wq = 0.5*wg

_cache = {}


def Nfun_of(f):
    f3 = sp.diff(f, x, 3); f4 = sp.diff(f, x, 4)
    F = sp.lambdify(x, f, 'numpy')
    F3 = sp.lambdify(x, f3, 'numpy')
    F3v = sp.lambdify(x, f3, 'numpy')
    v3g = F3(xg)  # array

    def N(xv):
        f3x = float(F3(xv))
        p3v = v3g + 0.0
        ig1 = (p3v - f3x)/(xv - xg)
        K2 = 1.0/(xv+xg) + 2.0*xg*(xv-xg)/(xv+xg)**3
        return float(np.dot(wq, ig1) + f3x*np.log(xv/(1.0-xv)) - np.dot(wq, p3v*K2))
    return F, N


def build(fs, m):
    n = len(fs)
    nodes = [0.5*(1-np.cos(np.pi*(i+1)/(m+1))) for i in range(m)]
    On = np.zeros((m, n)); Vn = np.zeros((m, n))
    gen = [Nfun_of(f) for f in fs]
    for i, xi in enumerate(nodes):
        for k, (F, N) in enumerate(gen):
            On[i, k] = N(xi)
            Vn[i, k] = float(F(xi))
    return nodes, On, Vn


def eigs(fs, m):
    nodes, On, Vn = build(fs, m)
    L = np.linalg.solve(Vn, On)
    ev = np.linalg.eigvals(L)
    return sorted([-e.real for e in ev if abs(e.imag) < 1e-6*max(1, abs(e)) and -e.real > 1e-6])


def poly_phi(N):
    fs = [3*x - 3*x**2 + x**3]
    for j in range(N+1):
        fs.append(x*(1-x)**3*sp.legendre(j, 2*x-1))
    return fs


print("=== polynomial-only (Legendre) ===", flush=True)
for N in [6, 10, 14, 18, 22]:
    e = eigs(poly_phi(N), N+2)
    print("N=%2d n=%2d  %s" % (N, N+2, " ".join("%.4f" % v for v in e[:3])), flush=True)

print("=== + log functions ===", flush=True)
for N in [6, 10, 14, 18]:
    fs = poly_phi(N) + [x*sp.log(x), x**2*sp.log(x), x**3*sp.log(x),
                        (1-x)**4*sp.log(1-x), (1-x)**5*sp.log(1-x)]
    e = eigs(fs, len(fs))
    print("N=%2d n=%2d  %s" % (N, len(fs), " ".join("%.4f" % v for v in e[:3])), flush=True)