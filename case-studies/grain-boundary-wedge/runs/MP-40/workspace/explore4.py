import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import polynomial as npoly
from scipy.linalg import eig

gl_x, gl_w = leggauss(600)
XI = 0.5*(gl_x+1.0)
WQ = 0.5*gl_w


def S_integral_poly(cg, x):
    gx = npoly.polyval(x, cg)
    n = len(cg)
    qc = np.zeros(n-1)
    for j in range(1, n):
        for m in range(j):
            qc[m] += cg[j]*x**(j-1-m)
    intq = sum(qc[m]/(m+1) for m in range(len(qc))) if len(qc) else 0.0
    S1 = -intq + gx*np.log(abs(x/(1.0-x)))
    gxi = npoly.polyval(XI, cg)
    Sreg = -1.0/(x+XI) - 2.0*XI*(x-XI)/(x+XI)**3
    return S1 + np.sum(WQ*gxi*Sreg)


def deriv_poly(c, k):
    c = np.array(c, float)
    for _ in range(k):
        c = np.polynomial.polynomial.polyder(c)
    return c


def cheb_basis_poly(k):
    Tk = np.polynomial.chebyshev.cheb2poly([0]*k + [1])
    p = np.array([0.0])
    for i, c in enumerate(Tk):
        term = np.array([c])
        for _ in range(i):
            term = np.convolve(term, np.array([-1.0, 2.0]))
        if len(term) > len(p):
            p = np.concatenate([p, np.zeros(len(term)-len(p))])
        p[:len(term)] += term
    p = np.convolve(p, np.array([0.0, 1.0]))
    p = np.convolve(p, np.array([1.0, -3.0, 3.0, -1.0]))
    return p


def opmatrix(K):
    B = [cheb_basis_poly(k) for k in range(K)]
    Nmat = np.zeros((len(XI), K))
    Bmat = np.zeros((len(XI), K))
    for k in range(K):
        c3 = deriv_poly(B[k], 3)
        Nmat[:, k] = np.array([S_integral_poly(c3, x) for x in XI])
        Bmat[:, k] = npoly.polyval(XI, B[k])
    G = (Bmat*WQ[:, None]).T @ Nmat
    H = (Bmat*WQ[:, None]).T @ Bmat
    return G, H, B, Bmat


K = 24
G, H, B, Bmat = opmatrix(K)
print("||G-G^T||/||G|| =", np.linalg.norm(G-G.T)/np.linalg.norm(G))
print("||H-H^T||/||H|| =", np.linalg.norm(H-H.T)/np.linalg.norm(H))

ev, evec = eig(G, -H)
order = np.argsort(ev.real)
xs = np.linspace(0, 1, 401)
for idx in order[:6]:
    lam = ev[idx]
    a = evec[:, idx]
    coeff = np.zeros(len(B[0]), complex)
    for k in range(K):
        if len(B[k]) > len(coeff):
            coeff = np.concatenate([coeff, np.zeros(len(B[k])-len(coeff))])
        c = B[k]
        coeff[:len(c)] += a[k]*c
    prof = npoly.polyval(xs, coeff)
    m = np.max(np.abs(prof)) or 1
    prof = prof/m
    print("eval=%.6f%+.6fj  prof at x=.1,.25,.5,.75,.9: %s" % (
        lam.real, lam.imag, np.round([prof[np.argmin(abs(xs-v))] for v in [0.1, 0.25, 0.5, 0.75, 0.9]], 4)))