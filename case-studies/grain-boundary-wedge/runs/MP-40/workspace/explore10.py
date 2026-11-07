import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig, null_space

NMU = 900
gl_x, gl_w = leggauss(NMU)
XI = 0.5*(gl_x+1.0)
WQ = 0.5*gl_w


def Dfun(x):
    return -np.log(abs(1-x)) + np.log(1+x) + 6*x/(1+x) - 2*x**2/(1+x)**2 - 4.0


def pint_S_g(g, x):
    gx = g(np.array([x]))[0]
    gxi = g(XI)
    diff = x - XI
    close = np.abs(diff) < 1e-9
    ratio = np.empty_like(diff)
    safe = ~close
    ratio[safe] = (gxi[safe]-gx)/diff[safe]
    if np.any(close):
        hh = 1e-5
        ratio[close] = -(g(np.array([x+hh]))[0]-g(np.array([x-hh]))[0])/(2*hh)
    integ = 8.0*x*XI**2*ratio/(x+XI)**3
    integ[close] = 0.0
    return np.sum(WQ*integ) + gx*Dfun(x)


def basis_ops(ps):
    """basis sigma_k = x^{p_k}, returns G,H, constr matrix C and M,B at nodes"""
    K = len(ps)
    Bq = np.zeros((len(XI), K)); Nq = np.zeros((len(XI), K))
    gmat = {}; gpmat = {}
    for i, p in enumerate(ps):
        c3 = p*(p-1)*(p-2)
        g = (lambda pp, cc: (lambda t: cc*t**(pp-3.0)))(p, c3)
        gp = (lambda pp, cc: (lambda t: cc*(pp-3.0)*t**(pp-4.0)))(p, c3)
        Bq[:, i] = XI**p
        Nq[:, i] = np.array([pint_S_g(g, x) for x in XI])
    G = (Bq*WQ[:, None]).T @ Nq
    H = (Bq*WQ[:, None]).T @ Bq
    C = np.zeros((2, K))
    for i, p in enumerate(ps):
        C[0, i] = p
        C[1, i] = p*(p-1)
    return G, H, C


def run(ps, label):
    G, H, C = basis_ops(ps)
    Z = null_space(C)
    G2 = Z.T @ G @ Z; H2 = Z.T @ H @ Z
    ev = eig(G2, -H2, right=False)
    re = np.sort(np.array([v.real for v in ev if abs(v.imag) < 1e-6*max(1, abs(v.real))
                           and v.real > 0]))
    print("%-22s K=%2d  lam1=%.4f  lam2=%.4f  lam3=%.1f" %
          (label, len(ps), re[0], re[1], re[2] if len(re) > 2 else -1))
    return re


print("Basis x^p, p = integers>=1 :")
for n in [6, 8, 10, 12]:
    ps = list(range(1, n+2))[:-1]
for n in [6, 8, 10, 12]:
    ps = list(np.arange(1.0, 1.0+n))
    run(ps, "int p")

print("\nBasis x^p, p = 3/2,2,5/2,3,... :")
for n in [6, 8, 10, 12, 14]:
    ps = list(1.5 + 0.5*np.arange(n))
    run(ps, "half p>=1.5")

print("\nBasis x^p, p = 1,3/2,2,5/2,... (mixed) :")
for n in [8, 10, 12, 14]:
    ps = list(1.0 + 0.5*np.arange(n))
    run(ps, "mixed p>=1")