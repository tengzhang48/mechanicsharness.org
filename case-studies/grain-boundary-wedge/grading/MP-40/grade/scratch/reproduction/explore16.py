import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig, null_space

NMU = 800
gl_x, gl_w = leggauss(NMU)
XI = 0.5*(gl_x+1.0); WQ = 0.5*gl_w


def Dfun(x):
    return -np.log(abs(1-x)) + np.log(1+x) + 6*x/(1+x) - 2*x**2/(1+x)**2 - 4.0


def pint_S_all(g, gp, xv):
    gx = g(xv); gxi = g(XI)
    diff = xv[:, None] - XI[None, :]
    close = np.abs(diff) < 1e-9
    num = gxi[None, :] - gx[:, None]
    ratio = np.where(close, -gp(xv)[:, None], num/np.where(close, 1.0, diff))
    integ = 8.0*xv[:, None]*XI[None, :]**2*ratio/(xv[:, None]+XI[None, :])**3
    integ = np.where(close, 0.0, integ)
    return integ @ WQ + gx*Dfun(xv)


def chebT_poly(k):
    return np.polynomial.chebyshev.cheb2poly([0]*k+[1.])


def to01(c):
    out = np.array([0.0])
    for i, ci in enumerate(c):
        term = np.array([ci])
        for _ in range(i):
            term = np.convolve(term, np.array([-1.0, 2.0]))
        if len(term) > len(out):
            out = np.concatenate([out, np.zeros(len(term)-len(out))])
        out[:len(term)] += term
    return out


def deriv(c, k):
    for _ in range(k):
        c = np.polynomial.polynomial.polyder(c)
    return c


def build(K):
    polys = [to01(chebT_poly(k)) for k in range(K+1)]
    C = np.zeros((3, K+1))
    for k, p in enumerate(polys):
        C[0, k] = np.polynomial.polynomial.polyval(0.0, p)
        C[1, k] = np.polynomial.polynomial.polyval(1.0, deriv(p, 1))
        C[2, k] = np.polynomial.polynomial.polyval(1.0, deriv(p, 2))
    Z = null_space(C)
    Bmat = np.zeros((len(XI), K+1)); Nmat = np.zeros((len(XI), K+1))
    for k, p in enumerate(polys):
        Bmat[:, k] = np.polynomial.polynomial.polyval(XI, p)
        c3 = deriv(p, 3)
        g = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, cc)))(c3)
        gp = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, deriv(cc, 1))))(c3)
        Nmat[:, k] = pint_S_all(g, gp, XI)
    Bn = Bmat @ Z; Nn = Nmat @ Z
    n = Bn.shape[1]; R = np.eye(n)
    for j in range(n):
        v = Bn[:, j].copy(); r = R[:, j].copy()
        for i in range(j):
            c = np.sum(WQ*Bn[:, i]*v); v -= c*Bn[:, i]; r -= c*R[:, i]
        R[:, j] = r/np.sqrt(np.sum(WQ*v*v))
    Q = Bn @ R; NQ = Nn @ R
    return (Q*WQ[:, None]).T @ NQ, (Q*WQ[:, None]).T @ Q, Q, NQ


def eigs(K):
    G, H, Q, NQ = build(K)
    ev, V = np.linalg.eig(np.linalg.solve(H, G))
    lam = [(-ev[k].real, k) for k in range(len(ev))
           if abs(ev[k].imag) < 1e-6*max(1, abs(ev[k])) and -ev[k].real > 0]
    lam.sort()
    return [t[0] for t in lam], ev, V, Q, NQ


for K in [20, 24, 28, 32, 36]:
    ls, *_ = eigs(K)
    print("K=%2d  " % K + "  ".join("%.4f" % v for v in ls[:5]), flush=True)