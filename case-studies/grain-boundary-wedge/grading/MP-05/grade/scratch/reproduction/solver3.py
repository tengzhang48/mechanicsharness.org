"""Solver 3: Galerkin eigen-solver using the sqrt-tip basis.

Basis: rho_n = sqrt(1-x) T_n(2x-1); phi_n = f[rho_n] (triple integral with BCs
f(0)=f'(1)=f''(1)=0 built in, phi_n''' = rho_n exactly).
sigma = sum c_n phi_n:  sigma_s = sum c_n K[rho_n]  (exact via core machinery).
Galerkin: <phi_k, sigma_s> = <phi_k, L sigma>  ->  M cdot = G c, M Gram, G proj.
Eigen of M^{-1} G  ->  mu = -lambda; filter by residual & N-convergence.
"""
import numpy as np
from core import kernel_matrix, f_matrix, cheb_basis, gauss_jacobi01


def galerkin_matrices(N, Mquad=280, Kleg=140, Mleg=240):
    xq, wq = np.polynomial.legendre.leggauss(Mleg)
    xq = (xq + 1) / 2
    wq = wq / 2
    A = kernel_matrix(xq, N, Mquad)            # (M, N+1): K[rho_n](xq)
    B = f_matrix(xq, N, Mquad, Kleg)           # (M, N+1): phi_n(xq)
    G = A.T @ (wq[:, None] * B)                # G[k,n] = <phi_k, K rho_n>
    M = B.T @ (wq[:, None] * B)                # M[k,n] = <phi_k, phi_n>
    return M, G


def galerkin_eigs(N, Mquad=280, Kleg=140, Mleg=240):
    M, G = galerkin_matrices(N, Mquad, Kleg, Mleg)
    from scipy.linalg import eig, solve
    X = solve(M, G)
    mu, Z = eig(X)
    return mu, Z, M, G, X
