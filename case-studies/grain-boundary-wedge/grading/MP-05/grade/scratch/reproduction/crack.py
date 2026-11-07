"""Steady state (t->inf): traction-free GB = surface crack under remote sigma0.
Solve K[rho_inf](x) = 1 on (0,1) with rho_inf = (1-x)^{-1/2} P(x) (crack tip).
K_inf = -pi*sqrt(2)*A, A = lim rho_tilde*sqrt(1-x) = P(1)."""
import numpy as np
import scipy.special as sp
from solver5 import gauss_jacobi01, cheb_T


def J0inv(x):
    d = np.sqrt(1.0 - np.asarray(x, dtype=float))
    return np.log(np.maximum(1.0-d, 1e-300)/(1.0+d))/d


def crack_ops(N, Mquad=280):
    xg, wg = gauss_jacobi01(Mquad, -0.5)
    Qg = cheb_T(xg, N)
    return xg, wg, Qg


def hinv_matrix(x, N, Mquad=280):
    """[x,n] = K[ (1-t)^{-1/2} T_n ](x)."""
    xg, wg, Qg = crack_ops(N, Mquad)
    Qc = cheb_T(x, N)
    A = Qc * J0inv(x)[:, None]
    num = Qg[None, :, :] - Qc[:, None, :]
    dxi = x[:, None] - xg[None, :]
    A = A + np.einsum('j,ijn->in', wg, num / dxi[:, :, None])
    k1 = -1.0/(x[:, None] + xg[None, :])
    k2 = -2.0*xg[None, :]*dxi/(x[:, None] + xg[None, :])**3
    A = A + np.einsum('j,ij,jn->in', wg, k1, Qg)
    A = A + np.einsum('j,ij,jn->in', wg, k2, Qg)
    return A


def solve_crack(N=32, Mquad=280):
    ii = np.arange(1, N + 2)
    xc = (1.0 - np.cos((2*ii-1)*np.pi/(2*(N+1))))/2.0
    A = hinv_matrix(xc, N, Mquad)
    a = np.linalg.solve(A, np.ones(N + 1))
    return a, xc, A


def stress_ahead(xa, a, N, Mquad=280):
    """sigma_total(x) for x>1: 1 - K[rho](x), no PV singularity."""
    xg, wg, Qg = crack_ops(N, Mquad)
    Qa = cheb_T(xa, N)
    rho_g = (1.0 - xg)**(-0.5) * (Qg @ a)
    rho_a = (1.0 - xa)**(-0.5) * (Qa @ a)
    ker = 1.0/(xa[:, None] - xg[None, :]) - 1.0/(xa[:, None] + xg[None, :]) \
          - 2.0*xg[None, :]*(xa[:, None] - xg[None, :])/(xa[:, None] + xg[None, :])**3
    return 1.0 - (ker * wg) @ rho_g


if __name__ == "__main__":
    for N in [16, 24, 32, 40, 48]:
        a, xc, A = solve_crack(N)
        res = np.abs(A @ a - 1).max()
        Acoef = float(np.sum(a))       # P(1) = sum a_n (T_n(1)=1)
        Knorm = -np.pi*np.sqrt(2)*Acoef
        print(f"N={N}: colloc residual={res:.2e}  A={Acoef:+.8f}  K_norm={Knorm:.8f}")
