import numpy as np
from numpy.polynomial.legendre import leggauss

NMU = 1200
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


def make_basis(M):
    # phi0 = x*(1 + w + w^2),  phi_j = x*w^j (j>=3),  w=1-x
    basis = []
    w = np.array([1.0, -1.0])
    one_plus = np.polyadd(np.array([1.0]), w)
    one_plus = np.polyadd(one_plus, np.convolve(w, w))
    basis.append(np.convolve(np.array([0.0, 1.0]), one_plus))
    for j in range(3, M+1):
        pw = np.array([1.0])
        for _ in range(j):
            pw = np.convolve(pw, w)
        basis.append(np.convolve(np.array([0.0, 1.0]), pw))
    return basis


def deriv(c, k):
    for _ in range(k):
        c = np.polynomial.polynomial.polyder(c)
    return c


def run(M):
    basis = make_basis(M)
    K = len(basis)
    Bmat = np.zeros((len(XI), K)); Nmat = np.zeros((len(XI), K))
    for k, p in enumerate(basis):
        Bmat[:, k] = np.polynomial.polynomial.polyval(XI, p)
        c3 = deriv(p, 3)
        g = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, cc)))(c3)
        gp = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, deriv(cc, 1))))(c3)
        Nmat[:, k] = pint_S_all(g, gp, XI)
    # orthonormalize
    n = K; R = np.eye(n)
    for j in range(n):
        v = Bmat[:, j].copy(); r = R[:, j].copy()
        for i in range(j):
            c = np.sum(WQ*Bmat[:, i]*v); v -= c*Bmat[:, i]; r -= c*R[:, i]
        R[:, j] = r/np.sqrt(np.sum(WQ*v*v))
    Q = Bmat @ R; NQ = Nmat @ R
    G = (Q*WQ[:, None]).T @ NQ; H = (Q*WQ[:, None]).T @ Q
    ev = np.linalg.eigvals(np.linalg.solve(H, G))
    lam = [(-e.real, e) for e in ev if abs(e.imag) < 1e-6*max(1, abs(e)) and -e.real > 0]
    lam.sort(key=lambda t: t[0])
    return [t[0] for t in lam], Q, NQ, ev


for M in [8, 12, 16, 20, 24, 28]:
    ls, Q, NQ, ev = run(M)
    print("M=%2d " % M + "  ".join("%.5f" % v for v in ls[:5]), flush=True)