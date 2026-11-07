"""Solver 8 (fixed): resolvent with enriched basis.

sigma_hat(x;p) = alpha0*f0(x) + sum a_n f_n(x), where
  rho_inf (crack density, exact from crack.py)      -> f0 = f[rho_inf]
  rho_n = (1-x)^{3/2} T_n                            -> f_n
f[rho] = x*I1 - (x^2/2)*I0 + (1/2) int_0^x (x-t)^2 rho(t) dt  (BCs built in).
Equation rows at collocation x_i:   p*f[rho](x_i) - K[rho](x_i) = 1
with K[rho_inf] = 1 (exact), K[rho_n] = h_n (exact machinery).
Poles of sigma_hat in p = -lambda_n (decay rates); residues = mode shapes.
"""
import numpy as np
from solver7 import RhoOps32
from solver5 import cheb_T, gauss_jacobi01
from crack import solve_crack


class Resolvent:
    def __init__(self, N, Mquad=280, Kleg=120, Ncrack=40):
        self.N = N
        self.ops = RhoOps32(N, Mquad)
        # crack component
        ac, _, _ = solve_crack(Ncrack)
        self.acrack = ac
        # moments for crack basis: weight (1-t)^{-1/2}
        xgm, wgm = gauss_jacobi01(Mquad, -0.5)
        Qgm = cheb_T(xgm, Ncrack)
        rho_inf_g = (1 - xgm)**(-0.5) * (Qgm @ ac)
        self.rho_inf_g = rho_inf_g
        self.I0_cr = wgm @ rho_inf_g
        self.I1_cr = wgm @ (xgm * rho_inf_g)
        self.Kcr = 1.0   # K[rho_inf](x) = 1 exactly
        # moments for (1-x)^{3/2} basis
        xg, wg = gauss_jacobi01(Mquad, 1.5)
        Qg = cheb_T(xg, N)
        self.I0 = wg @ ((1 - xg)[:, None] * Qg)
        self.I1 = wg @ (xg[:, None] * (1 - xg)[:, None] * Qg)
        t, w = np.polynomial.legendre.leggauss(Kleg)
        self.t = (t + 1) / 2; self.w = w / 2
        self.tg = xg; self.wg = wg; self.Qg = Qg
        # collocation
        ii = np.arange(1, N + 2)
        self.xc = (1.0 - np.cos((2 * ii - 1) * np.pi / (2 * (N + 1)))) / 2.0
        self.H = self.ops.h_pieces(self.xc, kmax=0)[0]     # K[rho_n]
        self.F = self._frows_basis(self.xc)                  # f_n(xc)
        self.F0 = self._frows_crack(self.xc)                 # f0(xc)

    def _frows_basis(self, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        F = np.zeros((x.size, self.N + 1))
        for i, xv in enumerate(x):
            F[i, :] = xv * self.I1 - 0.5 * xv * xv * self.I0
            if xv > 0:
                s = xv * self.t
                Qs = cheb_T(s, self.N)
                ws = self.w * (xv - s) ** 2 * (1 - s) ** 1.5
                F[i, :] += 0.5 * xv * (ws @ Qs)
        return F

    def _frows_crack(self, x):
        from numpy.polynomial.chebyshev import chebval
        x = np.atleast_1d(np.asarray(x, dtype=float))
        F = np.zeros(x.size)
        for i, xv in enumerate(x):
            F[i] = xv * self.I1_cr - 0.5 * xv * xv * self.I0_cr
            if xv > 0:
                s = xv * self.t
                rs = (1 - s)**(-0.5) * (cheb_T(s, len(self.acrack) - 1) @ self.acrack)
                F[i] += 0.5 * xv * (self.w @ ((xv - s) ** 2 * rs))
        return F

    def solve(self, p):
        A = np.hstack([(p * self.F0 - 1.0)[:, None], p * self.F - self.H])
        b = np.ones(self.N + 1)
        v = np.linalg.lstsq(A, b, rcond=None)[0]
        return v

    def sigmahat(self, p, xout=0.5):
        v = self.solve(p)
        fo = self._frows_basis(np.array([xout]))[0] @ v[1:]
        f0o = self._frows_crack(np.array([xout]))[0] * v[0]
        return fo + f0o, v

    def mode_from_residue(self, pstar, xout, eps=1e-6):
        """residue of sigma_hat at pole pstar: (p-pstar)*sigma_hat(p)."""
        pm = pstar * (1 - eps)
        v1, _ = self.solve(pm)
        sm, _ = self.sigmahat(pstar * (1 - eps), xout)
        sp, _ = self.sigmahat(pstar * (1 + eps), xout)
        return (pm - pstar * (1 + eps)) * 0 + (pm - pstar) * sm  # ~ residue
