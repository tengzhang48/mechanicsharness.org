"""Solver 9: projected evolution matrix in enriched basis.

State: rho = alpha*rho_inf + sum a_n (1-x)^{3/2} T_n(2x-1).
Rho_s = (K rho)''' exact rows; note (K rho_inf)''' = 0 exactly.
Weighted least-squares projection with SVD truncation -> reduced matrix E;
eigenvalues of E approximate -lambda_n; verify by residuals vs continuous
equation at off-grid points and by N-convergence.
"""
import numpy as np
from solver7 import RhoOps32
from solver5 import cheb_T, gauss_jacobi01
from crack import solve_crack


def build_E(N=24, Mquad=280, Mleg=200, Ncrack=40, rcond=1e-11):
    ops = RhoOps32(N, Mquad)
    ac, _, _ = solve_crack(Ncrack)
    # quadrature for projection (plain Gauss-Legendre on [0,1])
    t, w = np.polynomial.legendre.leggauss(Mleg)
    xq = (t + 1) / 2; wq = w / 2
    # basis values at xq
    B = (1 - xq)[:, None] ** 1.5 * cheb_T(xq, N)            # (M, N+1)
    rho_inf_q = (1 - xq) ** (-0.5) * (cheb_T(xq, Ncrack) @ ac)
    Bf = np.hstack([rho_inf_q[:, None], B])                  # (M, N+2)
    # T rows: (K rho_n)''' at xq (crack component -> 0)
    H3 = ops.h_pieces(xq, kmax=3)[3]                         # (M, N+1)
    T = np.hstack([np.zeros((Mleg, 1)), H3])                 # (M, N+2)
    # weighted LSQ: M0 cdot = proj of T c ; M0 = Bf^T W Bf
    M0 = Bf.T @ (wq[:, None] * Bf)
    rhs = Bf.T @ (wq[:, None] * T)
    U, s, Vt = np.linalg.svd(M0)
    keep = s > rcond * s[0]
    Uk = U[:, keep]; sk = s[keep]
    # E = M0^+ T-proj  restricted to kept subspace:
    E = (Uk / sk).T @ rhs                                    # (r, N+2)
    # map back: c = Uk y  -> evolution of y: ydot = E c = E Uk y
    Ey = E @ Uk
    return Ey, Uk, keep.sum(), ops, ac, xq, wq, Bf, T


def spectrum(Ey, nshow=8):
    mu, Z = np.linalg.eig(Ey)
    sel = [(-m.real, m, i) for i, m in enumerate(mu)
           if abs(m.imag) < 1e-8 * max(1, abs(m.real)) and m.real < 0]
    sel.sort(key=lambda t: -t[0] * 0 - t[0] * -1 * 0 + t[0] * 0)  # no-op
    sel.sort(key=lambda t: -t[0])
    return sel[:nshow], mu, Z


if __name__ == "__main__":
    import warnings; warnings.filterwarnings('ignore')
    for N in [16, 20, 24, 28, 32]:
        Ey, Uk, r, ops, ac, xq, wq, Bf, T = build_E(N)
        sel, mu, Z = spectrum(Ey)
        lams = [-m for _, m, _ in sel]
        print(f"N={N} rank={r}: slowest decay rates lambda: "
              + "  ".join(f"{l:.6f}" for l in sorted(lams)[:6]))
