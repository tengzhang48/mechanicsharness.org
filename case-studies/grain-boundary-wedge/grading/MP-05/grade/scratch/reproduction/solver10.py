"""Solver 10: resolvent in the crack-like basis rho_n = (1-x)^{-1/2} T_n.

All BCs on sigma_hat = f[rho_hat] are built in; K[rho_n] exact (crack
machinery).  Square system  p*f[rho](x_i) - K[rho](x_i) = 1; poles in p
= -lambda_n detected via smallest singular value.
"""
import numpy as np
from solver5 import cheb_T, gauss_jacobi01
from crack import hinv_matrix


class Res10:
    def __init__(self, N, Mquad=280, Kleg=120):
        self.N = N
        xg, wg = gauss_jacobi01(Mquad, -0.5)
        Qg = cheb_T(xg, N)
        self.I0 = wg @ Qg                       # int (1-t)^-1/2 T_n
        self.I1 = wg @ (xg[:, None] * Qg)
        t, w = np.polynomial.legendre.leggauss(Kleg)
        self.t = (t + 1) / 2; self.w = w / 2
        ii = np.arange(1, N + 2)
        self.xc = (1.0 - np.cos((2 * ii - 1) * np.pi / (2 * (N + 1)))) / 2.0
        self.H = hinv_matrix(self.xc, N, Mquad)
        self.F = self._frows(self.xc)

    def _frows(self, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        F = np.zeros((x.size, self.N + 1))
        for i, xv in enumerate(x):
            F[i, :] = xv * self.I1 - 0.5 * xv * xv * self.I0
            if xv > 0:
                s = xv * self.t
                Qs = cheb_T(s, self.N)
                ws = self.w * (xv - s) ** 2 * (1 - s) ** (-0.5)
                F[i, :] += 0.5 * xv * (ws @ Qs)
        return F

    def smin(self, p):
        A = p * self.F - self.H
        return np.linalg.svd(A, compute_uv=False)[-1]

    def solve(self, p):
        return np.linalg.solve(p * self.F - self.H, np.ones(self.N + 1))

    def sigmahat(self, p, xo):
        a = self.solve(p)
        return self._frows(np.array([xo]))[0] @ a


if __name__ == "__main__":
    import warnings; warnings.filterwarnings('ignore')
    from scipy.optimize import minimize_scalar
    for N in [16, 20, 24, 28]:
        R = Res10(N)
        # check p=0: should reproduce crack: K rho = -1... note eq: pF-H = 1 -> -Ha=1
        a0 = R.solve(0.0)
        ls = np.logspace(np.log10(0.05), np.log10(5000), 1200)
        v = np.array([R.smin(-l) for l in ls])
        dips = [i for i in range(1, len(ls) - 1) if v[i] < v[i-1] and v[i] < v[i+1]]
        found = []
        for i in dips:
            r = minimize_scalar(lambda l: np.log10(R.smin(-l)), bounds=(ls[i-1], ls[i+1]),
                                method='bounded', options={'xatol': 1e-11})
            if 10**r.fun < 0.05 * min(v[i-1], v[i+1]):
                found.append((r.x, 10**r.fun))
        print(f"N={N}: {len(found)} poles:", [(f"{l:.6f}", f"{s:.1e}") for l, s in found[:8]])
