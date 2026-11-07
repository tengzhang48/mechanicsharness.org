"""Solver 7: same as solver5 but with rho-basis rho_n=(1-x)^{3/2} T_n(2x-1)."""
import numpy as np
import scipy.special as sp
from math import factorial
from scipy.linalg import eig
from solver5 import gauss_jacobi01, cheb_T, cheb_dcoefs, clenshaw, qderiv_vals, Dker


def J0w(x):
    """PV int_0^1 (1-t)^{3/2}/(x-t) dt = 2/3 + 2d^2 + d^3 l, d=sqrt(1-x), l=ln((1-d)/(1+d))."""
    d = np.sqrt(1.0 - np.asarray(x, dtype=float))
    u = 1.0 - d**2
    ell = np.log(np.maximum(1.0 - d, 1e-300) / (1.0 + d))
    return 2.0/3.0 + 2.0*d**2 + d**3*ell


def J0w_series(x, kmax=3):
    from math import factorial as F
    u = 1.0 - np.asarray(x, dtype=float)
    # J0 = 2/3 + 2u + sum_m a_m u^m ; term-by-term: expand d^3 l = u^{3/2} * l... use
    # direct series in u: J0(u) = 2/3 + 2u - 2 sum_{m>=2} c_m u^m? Derive numerically instead:
    # use exact: d^3*l = -2 d^4 (1 + d^2/3 + d^4/5 + ...) = -2u^2 - (2/3)u^3 - (2/5)u^4 - ...
    M = 10
    out = []
    # J0 = 2/3 + 2u - 2*u^2 - (2/3)u^3 - (2/5)u^4 - ... - 2u^m/(2m-3) for m>=2
    for k in range(kmax + 1):
        val = np.zeros_like(u)
        for m in range(0, M + 1):
            if m == 0:
                a = 2.0/3.0
            elif m == 1:
                a = 2.0
            else:
                a = -2.0/(2.0*m - 3.0)
            if m >= k and a != 0.0:
                val += a * (F(m)/F(m-k)) * u**(m-k)
        out.append((-1.0)**k * val)
    return out


def J0w_pieces(x):
    x = np.asarray(x, dtype=float)
    d = np.sqrt(1.0 - x); u = x
    ell = np.log(np.maximum(1.0 - d, 1e-300)/(1.0 + d))
    J = 2.0/3.0 + 2.0*d**2 + d**3*ell
    # derivatives wrt x (hand-derived, verify by FD):
    # d/dx[d^3 l] etc; use symbolic results:
    t = d  # t = sqrt(u'), u'=1-x
    J1 = -4.0 + 3.0*t**2*ell/(2.0*t) * 0 + 0  # placeholder; replaced below
    # compute by explicit formulas:
    # J = 2/3 + 2 s + s^2 *? -- let s=1-x: J = 2/3 + 2s + s^{3/2} l(s), l=ln((1-ss)/(1+ss)), ss=sqrt(s)
    # dJ/dx = -dJ/ds ; dJ/ds = 2 + (3/2)s^{1/2} l + s^{3/2} l' ; l'(s) = dl/ds = -2/(s(1-s))
    s = 1.0 - x
    ss = np.sqrt(s)
    lp = -1.0/(np.sqrt(s)*(1.0-s))                      # dl/ds
    dJds = 2.0 + 1.5*ss*ell + s**1.5*lp
    J1 = -dJds
    lpp = 0.5*s**-1.5/(1.0-s) - ss**-1/(1.0-s)**2       # d2l/ds2
    d2 = 0.75*s**-0.5*ell + 1.5*ss*lp - 1.0/(1.0-s)**2
    J2 = d2
    lppp = -0.75*s**-2.5/(1.0-s) + s**-1.5/(1.0-s)**2 + 2.0*ss**-1/(1.0-s)**3
    d3ds3 = -0.375*s**-1.5*ell + 1.5*s**-0.5*lp + 1.5*ss*lpp - 2.0/(1.0-s)**3
    J3 = -d3ds3
    near = x > 0.9999
    if np.any(near):
        S = J0w_series(x[near], 3)
        J = np.where(near, S[0], J); J1 = np.where(near, S[1], J1)
        J2 = np.where(near, S[2], J2); J3 = np.where(near, S[3], J3)
    return J, J1, J2, J3


class RhoOps32:
    def __init__(self, N, Mquad=280):
        self.N = N
        self.xg, self.wg = gauss_jacobi01(Mquad, 1.5)
        self.Qg = cheb_T(self.xg, N)
        self.Pcoef = []
        for n in range(N + 1):
            m = max(n + 2, 4)
            xl = (1.0 - np.cos(np.arange(m + 1) * np.pi / m)) / 2.0
            DD = (cheb_T(xl, n)[:, n][:, None] - self.Qg[None, :, n]) / (xl[:, None] - self.xg[None, :])
            pv = self.wg @ DD.T
            V = cheb_T(xl, m)
            self.Pcoef.append(np.linalg.solve(V, pv))

    def Pderiv(self, x, kmax=3):
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
        J = J0w_pieces(x)
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
