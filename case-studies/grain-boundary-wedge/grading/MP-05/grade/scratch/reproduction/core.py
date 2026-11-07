"""Core machinery: singular-integral eigenproblem for GB traction relaxation.

Dimensionless evolution (x=z/h, s=t/t0):
  sigma_s = K[sigma'''],  K[g](x) = PV-int_0^1 S(x,xi) g(xi) dxi,
  S = 1/(x-xi) - 1/(x+xi) - 2 xi (x-xi)/(x+xi)^3,
  BC sigma(0)=0, sigma_x(1)=sigma_xx(1)=0, IC sigma=1.

Eigen-mode: sigma = f e^{-lam s}:  -lam f = K[f'''].
Represent rho = f''' = sqrt(1-x) p(x), p in span{T_0..T_N}(2x-1); then
  f(x) = x*int xi rho - (x^2/2)*int rho + (1/2) int_0^x (x-xi)^2 rho dxi
satisfies f(0)=0, f'(1)=f''(1)=0 identically and f'''=rho.
Collocate -lam f(x_i) = K rho(x_i) at Chebyshev-Gauss points -> generalized EVP.
"""
import numpy as np
import scipy.special as sp


def gauss_jacobi01(M, beta):
    ab = beta
    n = np.arange(M)
    an = (0.0 - ab * ab) / ((2 * n + ab) * (2 * n + ab + 2.0))
    nn = np.arange(1, M)
    bn = 4 * nn * (nn + ab) * nn * (nn + ab) / (((2 * nn + ab) ** 2) * ((2 * nn + ab) ** 2 - 1.0))
    b0 = 2.0 ** (ab + 1.0) / (ab + 1.0)
    J = np.diag(an) + np.diag(np.sqrt(bn), 1) + np.diag(np.sqrt(bn), -1)
    tv, V = np.linalg.eigh(J)
    w = b0 * V[0] ** 2
    return (tv + 1.0) / 2.0, w * 2.0 ** (-beta - 1.0)


def cheb_basis(x, N):
    y = 2.0 * np.asarray(x, dtype=float) - 1.0
    out = np.zeros((y.size, N + 1))
    out[:, 0] = 1.0
    if N >= 1:
        out[:, 1] = y
    for n in range(1, N):
        out[:, n + 1] = 2 * y * out[:, n] - out[:, n - 1]
    return out


def J0_sqrt(x):
    d = np.sqrt(1.0 - x)
    return 2.0 + d * np.log((1.0 - d) / (1.0 + d))


def J0_invsqrt(x):
    d = np.sqrt(1.0 - x)
    return np.log((1.0 - d) / (1.0 + d)) / d


def kernel_matrix(x_col, N, Mquad=280, weight='sqrt'):
    """A[i,n] = int_0^1 S(x_i,xi) rho_n(xi) dxi for rho_n = w(xi) T_n(2xi-1),
    w = sqrt(1-xi) or (1-xi)^{-1/2}. Exact for the PV part (closed form J0 +
    polynomial division), Gauss-Jacobi (Mquad) for the rest."""
    if weight == 'sqrt':
        nodes, wts = gauss_jacobi01(Mquad, 0.5)
        J0 = J0_sqrt(x_col)
    else:
        nodes, wts = gauss_jacobi01(Mquad, -0.5)
        J0 = J0_invsqrt(x_col)
    Qg = cheb_basis(nodes, N)          # (M, N+1)
    Qc = cheb_basis(x_col, N)          # (Nc, N+1)
    Xc = np.asarray(x_col)[:, None]
    Xn = nodes[None, :]
    dxi = Xc - Xn
    A = Qc * J0[:, None]
    num = Qg[None, :, :] - Qc[:, None, :]
    A = A + np.einsum('j,ijn->in', wts, num / dxi[:, :, None])
    k1 = -1.0 / (Xc + Xn)
    k2 = -2.0 * Xn * dxi / (Xc + Xn) ** 3
    A = A + np.einsum('j,ij,jn->in', wts, k1, Qg)
    A = A + np.einsum('j,ij,jn->in', wts, k2, Qg)
    return A


def f_matrix(x_col, N, Mquad=280, Kleg=140, x_col_extra=None):
    """B[i,n] = f_n(x_i) for rho_n = sqrt(1-xi) T_n(2xi-1) (with the BCs built in),
    plus optionally extra evaluation points appended as extra rows."""
    nodes, wts = gauss_jacobi01(Mquad, 0.5)
    Qg = cheb_basis(nodes, N)
    I0 = wts @ Qg                          # int_0^1 sqrt(1-xi) T_n
    I1 = wts @ (nodes[:, None] * Qg)       # int_0^1 xi sqrt(1-xi) T_n
    pts = list(np.atleast_1d(x_col))
    if x_col_extra is not None:
        pts += list(np.atleast_1d(x_col_extra))
    B = np.zeros((len(pts), N + 1))
    t, w = np.polynomial.legendre.leggauss(Kleg)
    t = (t + 1) / 2
    w = w / 2
    for i, x in enumerate(pts):
        B[i, :] = x * I1 - 0.5 * x * x * I0
        if x > 0:
            s = x * t
            Qs = cheb_basis(s, N)
            B[i, :] += 0.5 * x * (w @ (((x - s) ** 2 * np.sqrt(1 - s))[:, None] * Qs))
    return B


def colloc_points(Nc):
    i = np.arange(Nc)
    th = (2 * i + 1) * np.pi / (2.0 * Nc)
    return (1.0 - np.cos(th)) / 2.0


def solve_eig(N, Mquad=280, Kleg=140):
    from scipy.linalg import eig
    Nc = N + 1
    xc = colloc_points(Nc)
    A = kernel_matrix(xc, N, Mquad)
    B = f_matrix(xc, N, Mquad, Kleg)
    mu, V = eig(A, B)
    lam = -mu
    res = []
    for k in range(len(mu)):
        a = V[:, k]
        r = np.linalg.norm(A @ a - mu[k] * (B @ a)) / (np.linalg.norm(A @ a) + np.linalg.norm(B @ a) + 1e-300)
        res.append(r)
    return lam, V, np.array(res), A, B, xc
