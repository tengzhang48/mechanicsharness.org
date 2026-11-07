"""Solver 5: exact operator eigenproblem in rho-space (well-conditioned).

rho_n(x) = sqrt(1-x) T_n(2x-1); h_n(x) = K[rho_n](x) exactly:
  h_n = q_n J0 + P_n + Rat_n
  J0 = 2 + d ln((1-d)/(1+d)), d=sqrt(1-x)                      (PV part, closed)
  P_n = sum_{m=1}^n q_n^{(m)}(x) M_m(x)/m!,  M_m = int (1-t)^{1/2}(x-t)^{m-1} dt
  Rat_n = -int (1-t)^{1/2} q_n(t) D(x,t) dt,
  D = 1/(x+t)+2t(x-t)/(x+t)^3 and its x-derivatives in closed form.
Eigenproblem: rows (K rho)'''(x_i) = -lam rho(x_i) at collocation pts
plus BC rows (K rho)(0)=0, (K rho)'(1)=0, (K rho)''(1)=0.
"""
import numpy as np
import scipy.special as sp
from math import factorial
from scipy.linalg import eig


def gauss_jacobi01(M, beta):
    n = np.arange(M)
    an = (0.0 - beta * beta) / ((2 * n + beta) * (2 * n + beta + 2.0))
    nn = np.arange(1, M)
    bn = 4 * nn * (nn + beta) * nn * (nn + beta) / (((2 * nn + beta) ** 2) * ((2 * nn + beta) ** 2 - 1.0))
    b0 = 2.0 ** (beta + 1.0) / (beta + 1.0)
    Jm = np.diag(an) + np.diag(np.sqrt(bn), 1) + np.diag(np.sqrt(bn), -1)
    tv, V = np.linalg.eigh(Jm)
    w = b0 * V[0] ** 2
    return (tv + 1.0) / 2.0, w * 2.0 ** (-beta - 1.0)


def cheb_T(x, N):
    y = 2.0 * np.asarray(x) - 1.0
    out = np.zeros((y.size, N + 1))
    out[:, 0] = 1.0
    if N >= 1:
        out[:, 1] = y
    for n in range(1, N):
        out[:, n + 1] = 2 * y * out[:, n] - out[:, n - 1]
    return out


def cheb_dcoefs(N, kmax):
    """coeffs (in T_j(y)) of d^k/dy^k T_n for n=0..N -> dict k -> (N+1)x(N+k+1)."""
    base = np.zeros((N + 1, N + 1))
    for n in range(N + 1):
        base[n, n] = 1.0
    out = {0: base}
    for k in range(1, kmax + 1):
        prev = out[k - 1]
        c = np.zeros((N + 1, N + k + 1))
        for n in range(N + 1):
            d = np.zeros(prev.shape[1] + 2)
            for j in range(prev.shape[1] - 1, 0, -1):
                d[j - 1] = d[j + 1] + prev[n, j] * 2.0 * j
            d[0] *= 0.5
            nn = min(n - k + 1, c.shape[1])
            if nn > 0:
                c[n, :nn] = d[:nn]
        out[k] = c
    return out


def clenshaw(y, c):
    b1 = np.zeros_like(y); b2 = np.zeros_like(y)
    for j in range(len(c) - 1, -1, -1):
        b0 = 2 * y * b1 - b2 + c[j]
        b2 = b1; b1 = b0
    return b1 - y * b2


def qderiv_vals(x, N, kmax):
    """q_n^{(k)}(x), q_n = T_n(2x-1): returns [kmax+1, len(x), N+1]."""
    y = 2.0 * np.asarray(x) - 1.0
    dc = cheb_dcoefs(N, kmax)
    out = np.zeros((kmax + 1, len(x), N + 1))
    for k in range(kmax + 1):
        C = dc[k]
        for n in range(N + 1):
            out[k, :, n] = (2.0 ** k) * clenshaw(y, C[n])
    return out


def J0_series(x, kmax=3):
    """J0 derivatives near x=1 via convergent series in u=1-x.
    J0 = 2 - 2 sum_{m>=1} u^m/(2m-1); d/dx = -d/du."""
    from math import factorial as F
    u = 1.0 - np.asarray(x, dtype=float)
    M = 8
    out = []
    for k in range(kmax + 1):
        val = np.zeros_like(u)
        if k == 0:
            val += 2.0
        for j in range(M + 1):
            a = 0.0 if j == 0 else (-2.0 / (2.0 * j - 1.0))
            if j >= k and a != 0.0:
                val += a * (F(j) / F(j - k)) * u**(j - k)
        out.append((-1.0) ** k * val)
    return out


def J0_pieces(x):
    x = np.asarray(x, dtype=float)
    d = np.sqrt(1.0 - x)
    u = 1.0 - d**2
    ell = np.log(np.maximum((1.0 - d), 1e-300) / (1.0 + d))
    J = 2.0 + d * ell
    J1 = 1.0 / u - ell / (2.0 * d)
    J2 = -1.0 / (2.0 * d**2 * u) - ell / (4.0 * d**3) - 1.0 / u**2
    J3 = -3.0 / (4.0 * d**4 * u) + 1.0 / (2.0 * d**2 * u**2) - 3.0 * ell / (8.0 * d**5) + 2.0 / u**3
    # switch to series near x=1
    near = x > 0.9999
    if np.any(near):
        S = J0_series(x[near], 3)
        J = np.where(near, S[0], J)
        J1 = np.where(near, S[1], J1)
        J2 = np.where(near, S[2], J2)
        J3 = np.where(near, S[3], J3)
    return J, J1, J2, J3


def Mm_polys(N, x):
    """M_m and d/dx up to 3rd: returns [4, len(x), m-1] for m=1..N+1."""
    x = np.asarray(x, dtype=float)
    out = np.zeros((4, x.size, N + 1))
    for m in range(1, N + 2):
        k = m - 1
        for i in range(k + 1):
            c = float(sp.comb(k, i, exact=True)) * (-1.0)**i * sp.beta(i + 1, 1.5)
            p = k - i
            out[0, :, m - 1] += c * x**p
            if p >= 1: out[1, :, m - 1] += c * p * x**(p - 1)
            if p >= 2: out[2, :, m - 1] += c * p * (p - 1) * x**(p - 2)
            if p >= 3: out[3, :, m - 1] += c * p * (p - 1) * (p - 2) * x**(p - 3)
    return out


def Dker(x, xi, k):
    s = x[:, None] + xi[None, :]
    z = xi[None, :]
    if k == 0:
        return s**-1 + 2 * z * (x[:, None] - z) * s**-3
    if k == 1:
        return -s**-2 - 4 * z * s**-3 + 12 * z**2 * s**-4
    if k == 2:
        return 2 * s**-3 + 12 * z * s**-4 - 48 * z**2 * s**-5
    if k == 3:
        return -6 * s**-4 - 48 * z * s**-5 + 240 * z**2 * s**-6
    raise ValueError


class RhoOps:
    """Exact h_n^{(k)} = (K rho_n)^{(k)} for k=0..3."""

    def __init__(self, N, Mquad=280):
        self.N = N
        self.xg, self.wg = gauss_jacobi01(Mquad, 0.5)
        self.Qg = cheb_T(self.xg, N)
        # precompute P_n as polynomials (Chebyshev coeffs) via exact quadrature
        # of DD_n(x,t)=(q_n(x)-q_n(t))/(x-t) at Chebyshev-Lobatto points
        self.Pcoef = []
        for n in range(N + 1):
            m = max(n + 2, 4)
            xl = (1.0 - np.cos(np.arange(m + 1) * np.pi / m)) / 2.0
            DD = (cheb_T(xl, n)[:, n][:, None] - self.Qg[None, :, n]) / (xl[:, None] - self.xg[None, :])
            pv = self.wg @ DD.T                          # P_n(xl), length m+1
            V = cheb_T(xl, m)
            self.Pcoef.append(np.linalg.solve(V, pv))  # coeffs of T_k(2x-1)

    def Pderiv(self, x, kmax=3):
        """P_n^{(k)}(x): [kmax+1, len(x), N+1] via stable coefficient recurrence."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        y = 2 * x - 1
        out = np.zeros((kmax + 1, x.size, self.N + 1))
        for n in range(self.N + 1):
            c = self.Pcoef[n]
            for k in range(kmax + 1):
                cur = c
                for kk in range(1, k + 1):
                    dd = np.zeros(len(cur) + 1)
                    for j in range(len(cur) - 1, 0, -1):
                        dd[j - 1] = dd[j + 1] + 2.0 * j * cur[j]
                    dd[0] *= 0.5
                    cur = dd
                out[k, :, n] = (2.0 ** k) * clenshaw(y, cur)
        return out

    def h_pieces(self, x, kmax=3):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        N = self.N
        qd = qderiv_vals(x, N, kmax + 3)
        J = J0_pieces(x)
        Pd = self.Pderiv(x, kmax)
        out = np.zeros((kmax + 1, x.size, N + 1))
        for n in range(N + 1):
            for k in range(kmax + 1):
                for j in range(k + 1):
                    cf = float(sp.comb(k, j, exact=True))
                    out[k, :, n] += cf * qd[j, :, n] * J[k - j]
        out -= Pd
        for k in range(kmax + 1):
            Dk = Dker(x, self.xg, k)
            out[k] -= (Dk * self.wg[None, :]) @ self.Qg
        return out


def eig_GBeigen(N, Mquad=280, x0c=1e-7, x1c=None):
    """Assemble and solve: rows (K rho)'''(x_i) = -lam rho(x_i) + BC rows."""
    from scipy.linalg import eig as _eig
    ops = RhoOps(N, Mquad)
    # collocation points: Chebyshev-Gauss interior points of (0,1)
    Nc = N - 2
    ii = np.arange(1, Nc + 1)
    xc = (1.0 - np.cos((2 * ii - 1) * np.pi / (2 * Nc))) / 2.0
    hp = ops.h_pieces(xc, kmax=3)[3]          # [x, n] h_n'''
    rho = np.sqrt(1 - xc)[:, None] * cheb_T(xc, N)   # [x, n]
    # BC rows
    x0 = np.array([x0c])
    h0 = ops.h_pieces(x0, kmax=0)[0]          # [1, n]
    x1 = np.array([1.0 - 1e-12]) if x1c is None else np.array([x1c])
    h1 = ops.h_pieces(x1, kmax=2)
    A = np.zeros((N + 1, N + 1))
    B = np.zeros((N + 1, N + 1))
    A[:Nc, :] = hp
    B[:Nc, :] = rho
    A[Nc, :] = h0[0]
    A[Nc + 1, :] = h1[1, 0]
    A[Nc + 2, :] = h1[2, 0]
    mu, Z = _eig(A, B)
    return mu, Z, ops, xc


def eig_residual(a, lam, ops, xr):
    """continuous residual of eigenpair at points xr."""
    hp = ops.h_pieces(xr, kmax=3)[3]
    rho = np.sqrt(1 - xr)[:, None] * cheb_T(xr, ops.N)
    r = hp @ a + lam * (rho @ a)
    return np.abs(r).max() / max(np.abs(hp @ a).max(), 1e-300)
