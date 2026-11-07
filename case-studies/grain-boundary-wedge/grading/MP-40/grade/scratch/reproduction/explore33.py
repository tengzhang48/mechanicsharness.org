exec(open('<round-runs>/runs/MP-40/explore26.py').read().split('prev = {}')[0])

import numpy as np
import mpmath as mp
from numpy.polynomial.legendre import leggauss


def poly_deriv(c, k):
    c = list(c)
    for _ in range(k):
        c = [c[i]*i for i in range(1, len(c))] if len(c) > 1 else [0.0]
    return c


def peval(c, x):
    r = 0.0
    for co in c[::-1]:
        r = r*x + co
    return r


Q = 160
tg, wg = leggauss(Q)
xg = 0.5*(tg+1.0)
wq = 0.5*wg


def Nphi(x, c):
    """Stable evaluation of N[phi](x)=int_0^1 s(x,xi) phi'''(xi) dxi."""
    c3 = poly_deriv(c, 3); c4 = poly_deriv(c, 4)
    p3x = peval(c3, x)
    # integrand 1: (p3(xi)-p3x)/(x-xi), smooth (limit -p3'(x))
    num = np.array([peval(c3, xi)-p3x for xi in xg])  # use polyval vectorized
    den = x - xg
    scale = max(abs(p3x), 1e-300)
    ig1 = num/den
    # term2
    K2 = 1.0/(x+xg) + 2.0*xg*(x-xg)/(x+xg)**3
    ig2 = np.array([peval(c3, xi) for xi in xg])*K2
    return float(np.dot(wq, ig1) + p3x*np.log(x/(1.0-x)) - np.dot(wq, ig2))


def Nphi_fast(x, c):
    c3 = poly_deriv(c, 3)
    p3v = np.polynomial.polynomial.polyval(xg, c3)
    p3x = np.polynomial.polynomial.polyval(x, c3)
    ig1 = (p3v - p3x)/(x - xg)
    K2 = 1.0/(x+xg) + 2.0*xg*(x-xg)/(x+xg)**3
    ig2 = p3v*K2
    return float(np.dot(wq, ig1) + p3x*np.log(x/(1.0-x)) - np.dot(wq, ig2))


# cross-check against sympy Nval for monomials
print("cross-check N[x^p] quadrature vs Nval at x=0.4,0.7:", flush=True)
for p in [4, 5, 6, 7, 8]:
    c = [0.0]*p + [1.0]
    for x in [0.4, 0.7]:
        q = Nphi_fast(x, c)
        s = float(Nval(p, mp.mpf(x)))
        print("  p=%d x=%.1f  quad=%.10f  nval=%.10f  rel=%.2e" % (p, x, q, s, abs(q-s)/max(1e-12, abs(s))), flush=True)