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
    Bmat = np.zeros((len(XI), K+1))
    Nmat = np.zeros((len(XI), K+1))
    for k, p in enumerate(polys):
        Bmat[:, k] = np.polynomial.polynomial.polyval(XI, p)
        c3 = deriv(p, 3)
        g = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, cc)))(c3)
        gp = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, deriv(cc, 1))))(c3)
        Nmat[:, k] = pint_S_all(g, gp, XI)
    Bn = Bmat @ Z
    Nn = Nmat @ Z
    # orthonormalize Bn in weighted L2, tracking R : Q = Bn @ R
    n = Bn.shape[1]
    R = np.eye(n)
    for j in range(n):
        v = Bn[:, j].copy(); r = R[:, j].copy()
        for i in range(j):
            c = np.sum(WQ*Bn[:, i]*v)
            v -= c*Bn[:, i]; r -= c*R[:, i]
        nrm = np.sqrt(np.sum(WQ*v*v))
        R[:, j] = r/nrm
    Q = Bn @ R
    NQ = Nn @ R
    G = (Q*WQ[:, None]).T @ NQ
    H = (Q*WQ[:, None]).T @ Q
    return G, H, Q, NQ


print("Correct BCs: sigma(0)=0, sigma'(1)=sigma''(1)=0, sigma(1) free")
prev = None
for K in [8, 12, 16, 20, 24]:
    G, H, Q, NQ = build(K)
    print("   G sym resid K=%d: %.2e" % (K, np.linalg.norm(G-G.T)/np.linalg.norm(G)))
    ev = eig(G, -H, right=False)
    re = np.sort(np.array([v.real for v in ev if abs(v.imag) < 1e-6*max(1, abs(v.real)) and v.real > 0]))
    print("K=%2d  lam1=%.6f lam2=%.5f lam3=%.4f lam4=%.3f" % (K, re[0], re[1], re[2], re[3]))