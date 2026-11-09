#!/usr/bin/env python3
"""Independent cross-checks of solve.py results:
   1. grid/tolerance sensitivity of all seven numbers (n_th x2, tighter rtol)
   2. |gamma| = 0.01 radii by root-finding on the field (not closed-form inversion)
   3. face slip by direct numerical integration of gamma_rz along the face ray
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from solve import solve_angular, interp

TAU0, GAM0, N, J = 80.0, 0.002, 3.5, 0.120


def run(n_th, rtol):
    # patch tolerance by recomputing here (solve_angular uses fixed rtol; grid is the main knob)
    lam, th, f, fp, T, I_n = solve_angular(N, n_th=n_th)
    p = 1.0 - lam
    c = GAM0 / TAU0 ** N
    A = (J / (c * I_n)) ** (1.0 / (N + 1.0))

    def Tv(q): return interp(th, T, q)
    def fv(q): return interp(th, f, q)
    def fpv(q): return interp(th, fp, q)

    def w(r, q): return c * A ** N * r ** p / p * Tv(q) ** (N - 1) * fpv(q)
    def tau_rz(r, q): return A * r ** (-p) * fpv(q)
    def tau_tz(r, q): return -A * lam * r ** (-p) * fv(q)
    def gam_abs(r, q): return GAM0 * (A * r ** (-p) * Tv(q) / TAU0) ** N

    out = {}
    out["dw20"] = 2 * w(0.0200, np.pi)
    out["dw80"] = 2 * w(0.0800, np.pi)
    out["trz"] = tau_rz(0.0200, np.pi / 3)
    out["ttz"] = tau_tz(0.0200, np.pi / 3)
    # root-find |gamma| = 0.01 on each ray (independent of closed-form inversion)
    for name, q in (("r0", 0.0), ("r90", np.pi / 2), ("r180", np.pi)):
        out[name] = brentq(lambda r: gam_abs(r, q) - 0.0100, 1e-6, 1.0, xtol=1e-14)
    # slip by integrating gamma_rz along the face ray theta -> pi
    def gamma_rz_face(r):
        tau = A * r ** (-p) * fpv(np.pi)          # T(pi) = fp(pi) > 0
        return c * tau ** N                        # gamma_rz = c tau^{n-1} tau_rz, tau_rz=tau
    for name, r in (("dw20_int", 0.0200), ("dw80_int", 0.0800)):
        val, _ = quad(gamma_rz_face, 0.0, r, epsabs=1e-15, epsrel=1e-12, limit=200)
        out[name] = 2 * val
    return lam, I_n, A, out


lam1, I1, A1, o1 = run(100001, 3e-13)
lam2, I2, A2, o2 = run(400001, 3e-13)

print(f"{'quantity':10s} {'grid 1e5':>18s} {'grid 4e5':>18s} {'rel diff':>12s}")
print(f"{'lambda':10s} {lam1:18.12f} {lam2:18.12f} {abs(lam1-lam2)/lam2:12.2e}")
print(f"{'I_n':10s} {I1:18.12f} {I2:18.12f} {abs(I1-I2)/abs(I2):12.2e}")
print(f"{'A':10s} {A1:18.12f} {A2:18.12f} {abs(A1-A2)/abs(A2):12.2e}")
for k in o1:
    print(f"{k:10s} {o1[k]:18.10e} {o2[k]:18.10e} {abs(o1[k]-o2[k])/abs(o2[k]):12.2e}")

print("\nclosed-form vs integration / root-find consistency:")
print(f"Delta w(0.02): closed {o2['dw20']:.10e} vs int {o2['dw20_int']:.10e} "
      f"rel {abs(o2['dw20']-o2['dw20_int'])/o2['dw20_int']:.2e}")
print(f"Delta w(0.08): closed {o2['dw80']:.10e} vs int {o2['dw80_int']:.10e} "
      f"rel {abs(o2['dw80']-o2['dw80_int'])/o2['dw80_int']:.2e}")
