import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig, null_space

NMU = 1500
gl_x, gl_w = leggauss(NMU)
XI = 0.5*(gl_x+1.0)
WQ = 0.5*gl_w


def Dfun(x):
    return -np.log(abs(1-x)) + np.log(1+x) + 6*x/(1+x) - 2*x**2/(1+x)**2 - 4.0


def pint_S_all(g, gp, xv):
    """P int_0^1 S(x,xi) g(xi) dxi for all x in xv, vectorized (g,gp array-capable)."""
    gx = g(xv)                 # (nx,)
    gxi = g(XI)                # (nmu,)
    diff = xv[:, None] - XI[None, :]
    close = np.abs(diff) < 1e-9
    num = gxi[None, :] - gx[:, None]
    number = np.where(close, -gp(xv)[:, None], num)
    ratio = np.where(close, number, num/np.where(close, 1.0, diff))
    integ = 8.0*xv[:, None]*XI[None, :]**2*ratio/(xv[:, None]+XI[None, :])**3
    integ = np.where(close, 0.0, integ)
    return integ @ WQ + gx*Dfun(xv)


def chebT_poly(k):
    return np.polynomial.chebyshev.cheb2poly([0]*k+[1.])


def to01(c):
    # compose polynomial c(t) with t = 2x-1  -> polynomial in x
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


def build_basis(K):
    polys = [to01(chebT_poly(k)) for k in range(K+1)]
    C = np.zeros((3, K+1))
    for k, p in enumerate(polys):
        C[0, k] = np.polynomial.polynomial.polyval(0.0, p)
        C[1, k] = np.polynomial.polynomial.polyval(1.0, deriv(p, 1))
        C[2, k] = np.polynomial.polynomial.polyval(1.0, deriv(p, 2))
    Z = null_space(C)                      # (K+1) x (K-2)
    # orthonormalize in L2(GL) : columns of Bmat
    Bmat = np.zeros((len(XI), K+1))
    for k, p in enumerate(polys):
        Bmat[:, k] = np.polynomial.polynomial.polyval(XI, p)
    Bn = Bmat @ Z
    # Gram-Schmidt (weighted)
    w = WQ
    Q = np.zeros_like(Bn)
    for j in range(Bn.shape[1]):
        v = Bn[:, j].copy()
        for i in range(j):
            v -= np.sum(w*Q[:, i]*v)*Q[:, i]
        v /= np.sqrt(np.sum(w*v*v))
        Q[:, j] = v
    return polys, Z, Q


def ops(polys, Z, Q):
    Kp = len(polys)
    Nmat = np.zeros((len(XI), Kp))
    for k, p in enumerate(polys):
        c3 = deriv(p, 3)
        g = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, cc)))(c3)
        gp = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, deriv(cc, 1))))(c3)
        Nmat[:, k] = pint_S_all(g, gp, XI)
    Nn = Nmat @ Z
    G = (Q*WQ[:, None]).T @ Nn
    H = (Q*WQ[:, None]).T @ Q
    return G, H


for K in [8, 12, 16, 20, 24, 28]:
    polys, Z, Q = build_basis(K)
    G, H = ops(polys, Z, Q)
    ev = eig(G, -H, right=False)
    re = np.sort(np.array([v.real for v in ev if abs(v.imag) < 1e-6*max(1, abs(v.real)) and v.real > 0]))
    print("K=%2d  lam1=%.6f lam2=%.5f lam3=%.4f lam4=%.3f" % (K, re[0], re[1], re[2], re[3]))