import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig, cholesky, solve_triangular

NMU = 1600
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


def dcoef(c, kk):
    for _ in range(kk):
        c = np.polynomial.polynomial.polyder(c)
    return c


def opmatrix(K):
    Bq = np.zeros((len(XI), K))
    Nq = np.zeros((len(XI), K))
    for k in range(K):
        c = cheb_basis_poly(k)
        c3 = dcoef(c, 3)
        Bq[:, k] = np.polynomial.polynomial.polyval(XI, c)
        g = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, cc)))(c3)
        Nq[:, k] = np.array([pint_S_g(g, x) for x in XI])
    G = (Bq*WQ[:, None]).T @ Nq
    H = (Bq*WQ[:, None]).T @ Bq
    return G, H


def solve(G, H):
    H = 0.5*(H+H.T)
    c = cholesky(H, lower=True)
    Gi = solve_triangular(c, G, lower=True)
    M = solve_triangular(c, Gi.T, lower=True).T
    ev = np.linalg.eigvals(M)
    return np.sort(ev.real)


for K in [8, 12, 16, 20]:
    G, H = opmatrix(K)
    ev = solve(G, H)
    print("K=%2d  first 6 eigs:" % K, np.round(ev[:6], 5))
print()
# high precision check that G is (nearly) symmetrizable: eigenvalues should be real
G, H = opmatrix(16)
c = cholesky(0.5*(H+H.T), lower=True)
Gi = solve_triangular(c, G, lower=True)
M = solve_triangular(c, Gi.T, lower=True).T
ev = np.linalg.eigvals(M)
print("max imag of eigenvalues (K=16):", np.max(np.abs(ev.imag)))
print("M symmetry resid (K=16):", np.linalg.norm(M-M.T)/np.linalg.norm(M))