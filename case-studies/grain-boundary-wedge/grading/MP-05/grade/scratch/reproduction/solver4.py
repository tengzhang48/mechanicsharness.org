"""Solver 4: robust time-evolution in the sqrt-tip Galerkin basis.

sigma = sum c_n phi_n (phi_n'''=rho_n=sqrt(1-x)T_n, BCs built in).
Galerkin ODE:  M cdot = G c,  M=<phi_k,phi_n>, G=<phi_k, K rho_n> (exact).
Regularize M by SVD truncation; integrate with expm; sigma(x,s)=sum c_n(s)phi_n(x).
Initial sigma0=1 projected in L^2.
"""
import numpy as np
from core import kernel_matrix, f_matrix
from scipy.linalg import expm


def build(N, Mquad=280, Kleg=140, Mleg=240, rcond=1e-12):
    xq, wq = np.polynomial.legendre.leggauss(Mleg)
    xq = (xq + 1) / 2
    wq = wq / 2
    A = kernel_matrix(xq, N, Mquad)
    B = f_matrix(xq, N, Mquad, Kleg)
    G = A.T @ (wq[:, None] * B)
    M = B.T @ (wq[:, None] * B)
    rhs = wq @ B                        # <phi_k, 1>
    Mp = np.linalg.pinv(M, rcond=rcond, hermitian=True)
    E = Mp @ G                          # cdot = E c  (regularized)
    c0 = Mp @ rhs
    return E, c0, B, xq


def evolve(N, ss, xout=0.5, rcond=1e-12, Mleg=240):
    E, c0, B, xq = build(N, rcond=rcond, Mleg=Mleg)
    bo = f_matrix(np.array([xout]), N)[0]
    vals = []
    c = c0
    for s in ss:
        c = expm(E * s) @ c0
        vals.append(bo @ c)
    return np.array(vals), E, c0, bo


def spectrum(E):
    mu, Z = np.linalg.eig(E)
    return mu, Z
