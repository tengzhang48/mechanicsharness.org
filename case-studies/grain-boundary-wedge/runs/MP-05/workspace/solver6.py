"""Solver 6: eigenvalues via least-squares rank loss (pollution-proof).

For trial lambda, build the overdetermined system (3Nc rows at Chebyshev-Gauss
points + 3 weighted BC rows):   [ (K rho_n)'''(x_i) + lam rho_n(x_i) ] a = 0,
BC: (K rho)(0) = 0, (K rho)'(1)=0, (K rho)''(1)=0.
Genuine eigenvalues = global minima of smallest singular value (near machine
zero). Eigenvector = corresponding right singular vector -> mode shape.
"""
import numpy as np
from solver5 import RhoOps, cheb_T


def build_Ml(ops, lam, Nc_mult=3, x0c=1e-7):
    N = ops.N
    Nc = Nc_mult * (N + 1)
    ii = np.arange(1, Nc + 1)
    xc = (1.0 - np.cos((2 * ii - 1) * np.pi / (2 * Nc))) / 2.0
    hp = ops.h_pieces(xc, kmax=3)[3]
    rho = np.sqrt(1 - xc)[:, None] * cheb_T(xc, N)
    M = np.zeros((Nc + 3, N + 1))
    M[:Nc, :] = hp + lam * rho
    h0 = ops.h_pieces(np.array([x0c]), kmax=0)[0]
    h1 = ops.h_pieces(np.array([1.0]), kmax=2)
    M[Nc, :] = h0[0]
    M[Nc + 1, :] = h1[1, 0]
    M[Nc + 2, :] = h1[2, 0]
    return M


def smin(ops, lam, Nc_mult=3):
    M = build_Ml(ops, lam, Nc_mult)
    s = np.linalg.svd(M, compute_uv=False)
    return s[-1], M


def refine(ops, lam0, Nc_mult=3, tol=1e-12, itmax=60):
    """Golden-section-ish local refinement of sigma_min(lam)."""
    from scipy.optimize import minimize_scalar
    def f(l):
        return np.log10(max(smin(ops, l, Nc_mult)[0], 1e-300))
    r = minimize_scalar(f, bracket=None, bounds=(lam0*0.7, lam0*1.4),
                        method='bounded', options={'xatol': 1e-10})
    return r.x, 10**r.fun
