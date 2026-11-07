import sympy as sp
import numpy as np
import mpmath

x, t = sp.symbols('x t', positive=True)

# D(x) = P int_0^1 8 x t^2 / ((x-t)(x+t)^3) dt  (the singular coefficient)
def clean(e):
    e = sp.expand(e)
    e = e.subs(sp.I*sp.pi, 0)
    e = sp.expand_log(e, force=True)
    return sp.simplify(sp.expand(e))


Dexpr = clean(sp.integrate(8*x*t**2/((x+t)**3*(x-t)), (t, 0, 1)))
print("D(x) =", Dexpr)

D = 18
# N[x^p](x) = p!/(p-3)! * [ int_0^1 8x t^2/(x+t)^3 * (t^{p-3}-x^{p-3})/(x-t) dt + x^{p-3} D(x) ]
Nfun = {}
for p in range(3, D+1):
    c = p*(p-1)*(p-2)
    reg = 8*x*t**2/(x+t)**3 * (t**(p-3)-x**(p-3))/(x-t)
    Ireg = clean(sp.integrate(reg, (t, 0, 1)))
    Np = c*(Ireg + x**(p-3)*Dexpr)
    Nfun[p] = clean(Np)
# N[x^p]=0 for p<3
for p in range(0, 3):
    Nfun[p] = sp.Integer(0)

print("N[x^3] =", Nfun[3])
# lambdify
Nnum = {p: sp.lambdify(x, Nfun[p], 'mpmath') for p in range(0, D+1)}


def evalN(p, v):
    return float(Nnum[p](mpmath.mpf(str(v))))


def make(D, m):
    nodes = 0.5*(1-np.cos(np.pi*np.arange(1, m+1)/(m+1)))  # interior Chebyshev
    Cb = np.zeros((3, D+1))
    for p in range(D+1):
        Cb[0, p] = 0.0**p if p == 0 else 0.0
        Cb[1, p] = p*1.0**(p-1) if p >= 1 else 0.0
        Cb[2, p] = p*(p-1)*1.0**(p-2) if p >= 2 else 0.0
    from scipy.linalg import null_space
    Z = null_space(Cb)
    V = np.zeros((m, D+1)); O = np.zeros((m, D+1))
    for p in range(D+1):
        V[:, p] = nodes**p
        for i in range(m):
            O[i, p] = evalN(p, nodes[i])
    VZ = V @ Z; OZ = O @ Z
    ev = np.linalg.eigvals(np.linalg.solve(VZ, OZ))
    lam = sorted([-e.real for e in ev if abs(e.imag) < 1e-6*max(1, abs(e)) and -e.real > 0])
    return lam, nodes, VZ, OZ


for D in [8, 10, 12, 14, 16, 18]:
    m = D+1-3
    lam, nodes, VZ, OZ = make(D, m)
    print("D=%2d m=%2d  " % (D, m) + "  ".join("%.4f" % v for v in lam[:5]))