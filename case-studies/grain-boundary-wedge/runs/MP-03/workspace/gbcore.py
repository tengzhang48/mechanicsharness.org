"""Core utilities: Chebyshev grids, differentiation, kernel quadrature (PV-aware)."""
import numpy as np
from numpy.polynomial.legendre import leggauss


def cheb_lobatto_01(N):
    """Chebyshev-Lobatto nodes on [0,1]: x[0]=0, x[N]=1."""
    k = np.arange(N + 1)
    return 0.5 * (1.0 - np.cos(np.pi * k / N))


def diff_mat(x):
    """General-nodes differentiation matrix via barycentric weights (exact for polynomials)."""
    n = len(x)
    w = (-1.0) ** np.arange(n)
    w[0] *= 0.5
    w[-1] *= 0.5
    D = np.zeros((n, n))
    for i in range(n):
        dx = x[i] - x
        D[i] = (w / w[i]) / dx
        D[i, i] = 0.0
        D[i, i] = -np.sum(D[i])
    return D


def bary_weights(x):
    """Barycentric weights for Chebyshev-Lobatto grid on [0,1]."""
    n = len(x) - 1
    w = (-1.0) ** np.arange(n + 1)
    w[0] *= 0.5
    w[-1] *= 0.5
    return w


def lagrange_eval(x, w, xi):
    """Matrix L[q, j] = L_j(xi_q) of Lagrange cardinal polynomials on grid x with weights w."""
    xi = np.atleast_1d(xi)
    diff = xi[:, None] - x[None, :]
    # handle exact node hits
    hit = diff == 0.0
    has_hit = hit.any(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        C = w[None, :] / diff
    L = C / np.sum(C, axis=1, keepdims=True)
    for q in np.where(has_hit)[0]:
        j = np.where(hit[q])[0][0]
        L[q] = 0.0
        L[q, j] = 1.0
    return L


def Rfun(x, xi):
    """Regular (non-PV) part of the kernel: Shat - 1/(x-xi)."""
    x = np.asarray(x, dtype=float)
    xi = np.asarray(xi, dtype=float)
    return -1.0 / (x + xi) - 2.0 * xi * (x - xi) / (x + xi) ** 3


def build_Q(xc, xgrid, M=None):
    """Matrix Q[i,j] = A[L_j](xc_i) = PV-int_0^1 Shat(xc_i,xi) L_j(xi) dxi.

    xc: collocation points (all in (0,1]); xgrid: grid defining cardinal functions L_j.
    Uses PV splitting: PV-int g(xi)/(x-xi) dxi = int [g(xi)-g(x)]/(x-xi) dxi + g(x) ln(x/(1-x)).
    """
    n = len(xgrid)
    if M is None:
        M = max(200, 8 * n)
    xi, wi = leggauss(M)
    xi = 0.5 * (xi + 1.0)
    wi = 0.5 * wi
    w = bary_weights(xgrid)
    L = lagrange_eval(xgrid, w, xi)  # (M, n)
    nc = len(xc)
    Q = np.zeros((nc, n))
    # delta_ij = L_j(xc_i) via cardinal property (evaluated below)
    # delta_ij = L_j(xc_i): cardinal property -> L(xc) if xc on grid else evaluate
    Lex = lagrange_eval(xgrid, w, xc)  # (nc, n)
    for i, x in enumerate(xc):
        d = Lex[i]  # L_j(x)
        t1 = (wi[:, None] * (L - d[None, :]) / (x - xi)[:, None]).sum(axis=0)
        t2 = d * np.log(x / (1.0 - x))
        t3 = (wi * Rfun(x, xi))[:, None] * L
        t3 = t3.sum(axis=0)
        Q[i] = t1 + t2 + t3
    return Q


def A_exact_const(x):
    """Exact A[1](x) = int_0^1 Shat(x,xi) dxi for x in (0,1)."""
    return (np.log(x**2 / (1.0 - x**2)) + 6 * x / (1 + x) + 2 * np.log(1 + x)
            - 2 * x**2 / (1 + x) ** 2 - 4.0 - 2 * np.log(x))
