import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import polynomial as npoly
from scipy.linalg import eig

# Gauss-Legendre nodes on [0,1] for the regular part
gl_x, gl_w = leggauss(400)
XI = 0.5*(gl_x+1.0)
W  = 0.5*gl_w


def S_integral_poly(cg, x):
    # P int_0^1 S(x,xi) g(xi) dxi for polynomial g with coeffs cg (increasing)
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
    return S1 + np.sum(W*gxi*Sreg)


def deriv_poly(c, k):
    c = np.array(c, float)
    for _ in range(k):
        c = np.polynomial.polynomial.polyder(c)
    return c


def basis_poly(k):
    # phi_k = x^(k+1) (1-x)^3
    c = np.array([1.0, -3.0, 3.0, -1.0])   # (1-x)^3
    c = np.concatenate([np.zeros(k+1), c])
    return c


def build(K, N):
    # collocation nodes = Chebyshev-Lobatto interior points on (0,1)
    x = 0.5*(1.0+np.cos(np.pi*np.arange(N+1)/N))
    xi_nodes = x[1:-1]  # interior only (exclude 0 and 1)
    Phi = np.zeros((len(xi_nodes), K))
    A = np.zeros((len(xi_nodes), K))
    for k in range(K):
        c = basis_poly(k)
        Phi[:, k] = npoly.polyval(xi_nodes, c)
        c3 = deriv_poly(c, 3)
        for i, xx in enumerate(xi_nodes):
            A[i, k] = S_integral_poly(c3, xx)
    return xi_nodes, Phi, A


for K in [12, 20, 30, 40]:
    N = 200
    xi, Phi, A = build(K, N)
    # Galerkin projection with test = basis samples (discrete)
    G = Phi.T @ A
    H = Phi.T @ Phi
    ev, evec = eig(G, -H)
    ev = np.real_if_close(ev)
    pos = np.sort(np.array([v.real for v in ev if abs(v.imag) < 1e-6*max(1, abs(v.real)) and v.real > 0]))
    print("K=%d  first 6 positive rates:" % K, np.round(pos[:6], 5))