import sympy as sp
import numpy as np
import mpmath as mp
from scipy.linalg import null_space

mp.mp.dps = 60
x, t = sp.symbols('x t', positive=True)


def clean(e):
    e = sp.expand(e); e = e.subs(sp.I*sp.pi, 0)
    e = sp.expand_log(e, force=True)
    return sp.simplify(sp.expand(e))


Dexpr = clean(sp.integrate(8*x*t**2/((x+t)**3*(x-t)), (t, 0, 1)))
DMAX = 44
Nfun = {}
for p in range(3, DMAX+1):
    c = p*(p-1)*(p-2)
    reg = 8*x*t**2/(x+t)**3*(t**(p-3)-x**(p-3))/(x-t)
    Ireg = clean(sp.integrate(reg, (t, 0, 1)))
    Nfun[p] = clean(c*(Ireg + x**(p-3)*Dexpr))
for p in range(0, 3):
    Nfun[p] = sp.Integer(0)
Nnum = {p: sp.lambdify(x, Nfun[p], 'mpmath') for p in range(DMAX+1)}
print("symbolic operator ready", flush=True)


def Nval(p, v):
    return Nnum[p](mp.mpf(v))


# Chebyshev T_k(2x-1) monomial coefficients (integer)
Tmat = np.zeros((DMAX+1, DMAX+1))
for k in range(DMAX+1):
    # T_k(y) with y=2x-1 ; coefficients in x via numpy chebyshev->poly
    cy = np.polynomial.chebyshev.cheb2poly([0]*k+[1.])
    # substitute y=2x-1
    out = np.array([0.0])
    for i, ci in enumerate(cy):
        term = np.array([ci])
        for _ in range(i):
            term = np.convolve(term, np.array([-1.0, 2.0]))
        if len(term) > len(out):
            out = np.concatenate([out, np.zeros(len(term)-len(out))])
        out[:len(term)] += term
    Tmat[:len(out), k] = out


def run(D):
    m = D-2
    nodes = [mp.mpf(0.5)*(1-mp.cos(mp.pi*(i+1)/(m+1))) for i in range(m)]
    O = mp.matrix(m, D+1); V = mp.matrix(m, D+1)
    for i in range(m):
        xi = nodes[i]
        for p in range(D+1):
            O[i, p] = Nval(p, xi)
            V[i, p] = xi**p
    T = mp.matrix([[mp.mpf(Tmat[r, k]) for k in range(D+1)] for r in range(D+1)])
    Och = O*T; Vch = V*T
    On = np.array([[float(Och[i, k]) for k in range(D+1)] for i in range(m)])
    Vn = np.array([[float(Vch[i, k]) for k in range(D+1)] for i in range(m)])
    # BC in chebyshev basis: sigma(0)=0, sigma'(1)=0, sigma''(1)=0  (x=1 -> y=1, T_k(1)=1, T_k'(1)=k^2, T_k''(1)=...)
    Cb = np.zeros((3, D+1))
    for k in range(D+1):
        Cb[0, k] = np.polynomial.chebyshev.chebval(-1.0, [0]*k+[1.])   # T_k(-1)=(-1)^k
        Cb[1, k] = 2.0* (k*k) if k >= 1 else 0.0   # d/dx T_k(2x-1)=2 T_k'(y), T_k'(1)=k^2
    # second derivative: d2/dx2 = 4 T_k''(y); T_k''(1) = (k^2(k^2-1))/3
    for k in range(D+1):
        Cb[2, k] = 4.0*(k*k*(k*k-1))/3.0
    Z = null_space(Cb)
    A = On @ Z; B = Vn @ Z
    ev = np.linalg.eigvals(np.linalg.solve(B, A))
    lam = sorted([-e.real for e in ev if abs(e.imag) < 1e-5*max(1, abs(e)) and -e.real > 1e-6])
    return lam, nodes, Z, On, Vn


for D in [8, 10, 12, 14, 16, 18, 20]:
    lam, nodes, Z, On, Vn = run(D)
    print("D=%2d  " % D + "  ".join("%.4f" % v for v in lam[:5]), flush=True)