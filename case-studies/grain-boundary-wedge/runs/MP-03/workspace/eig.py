"""Eigenanalysis of -lambda f = A[f'''], f(0)=0, f'(1)=f''(1)=0 (s-space collocation pencil)."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from gbcore import cheb_lobatto_01, diff_mat, bary_weights, lagrange_eval, Rfun


def build_Q_general(xc, xgrid, M=None):
    """Q[i,j] = A[L_j](xc_i) for xc in (0,1); PV split with log term."""
    n = len(xgrid)
    if M is None:
        M = max(400, 8 * n)
    xi, wi = leggauss(M)
    xi = 0.5 * (xi + 1.0)
    wi = 0.5 * wi
    w = bary_weights(xgrid)
    L = lagrange_eval(xgrid, w, xi)      # (M, n)
    Lex = lagrange_eval(xgrid, w, xc)    # (nc, n)
    Q = np.zeros((len(xc), n))
    for i, x in enumerate(xc):
        d = Lex[i]
        t1 = (wi[:, None] * (L - d[None, :]) / (x - xi)[:, None]).sum(axis=0)
        t2 = d * np.log(x / (1.0 - x))
        t3 = ((wi * Rfun(x, xi))[:, None] * L).sum(axis=0)
        Q[i] = t1 + t2 + t3
    return Q


def eigen_s(N=60, M=None, nkeep=6, verbose=True):
    xg = cheb_lobatto_01(N)          # x[0]=0 ... x[N]=1
    D = diff_mat(xg)
    D2 = D @ D
    D3 = D2 @ D
    Q = build_Q_general(xg[1:N - 1], xg, M=M)   # rows at x_1..x_{N-2}
    G = -(Q @ D3)                    # -A[f'''] at those rows -> G f = lambda f
    Lmat = np.zeros((N + 1, N + 1))
    Bmat = np.zeros((N + 1, N + 1))
    # PDE rows: i = 1..N-2 -> matrix rows 1..N-2
    Lmat[1:N - 1, :] = G
    Bmat[1:N - 1, 1:N - 1] = np.eye(N - 2)
    # BC rows
    Lmat[0, 0] = 1.0                       # f(0)=0
    Lmat[N - 1, :] = D2[N, :]              # f''(1)=0 (placed in row N-1)
    Lmat[N, :] = D[N, :]                   # f'(1)=0 (row N)
    import scipy.linalg as sla
    lam, V = sla.eig(Lmat, Bmat, right=True)
    lam = lam[np.isfinite(lam)]
    sel = (np.abs(lam.imag) < 1e-8 * np.maximum(np.abs(lam.real), 1)) & (lam.real > 1e-6)
    lams = np.sort(lam[sel].real)
    if verbose:
        print(f"N={N}: smallest real positive eigenvalues: {lams[:nkeep]}")
    return xg, lam, V, lams


if __name__ == "__main__":
    for N in [24, 32, 40, 48, 60, 72]:
        eigen_s(N=N)
