"""Independent recomputation of the V2 key (bilinear membrane, inverse for T), audit 2026-11-07.

Written from V2/TASK_PROBLEM.md only. Differs from the builder's solve_reference.py route A (closed
formulas nu_t, J = K^2/T, Kolosov displacement) and from the case page's hrr_check.py (same closed
formulas) in the steps at risk:
  * constitutive reduction: the leading compliance is obtained NUMERICALLY as lim eps(L*sigma)/L of
    the complete stated law (mpmath, L = 1e25) on the 9 basis stresses; no nu_t formula is used;
  * J normalization: J is evaluated NUMERICALLY on square contours, W = sigma:eps/2 from the full 3D
    tensors, u_{y,x} = eps_xy + omega with the rotation omega integrated along the contour from the
    compatibility equations (omega = 0 on the ligament); no J = K^2/T;
  * stresses: Westergaard function Z = K/sqrt(2 pi z) in Cartesian form;
  * C_y: by integrating eps_yy along x = 0 from the tip; no displacement formula;
  * T: brentq on the numerically computed C_zz(T); screens by root finding on sigma_e and on the
    eigenvalues of the 3D strain tensor.
Run: python -I recompute_V2.py   (writes recompute_V2_results.json next to this file)
"""
import json
import os
import numpy as np
import mpmath as mp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
E, nu, sy, Jgiven, Czz_given = 70000.0, 0.30, 10.0, 0.800, -0.005
mp.mp.dps = 60


def full_law_mp(sig, T):
    """Complete stated total-strain law (3x3 mpmath matrices)."""
    sk = sig[0, 0] + sig[1, 1] + sig[2, 2]
    I = mp.eye(3)
    s = sig - sk / 3 * I
    se = mp.sqrt(mp.mpf(3) / 2 * sum(s[i, j] ** 2 for i in range(3) for j in range(3)))
    fac = max(1 - sy / se, 0) if se > 0 else 0
    return (1 + nu) / E * sig - nu / E * sk * I + mp.mpf(3) / 2 * (1 / mp.mpf(T) - 1 / mp.mpf(E)) * fac * s


def full_law(sig, T):
    """Float version of the complete law (used for a mutant only)."""
    sk = np.trace(sig); s = sig - sk / 3 * np.eye(3)
    se = np.sqrt(1.5 * np.sum(s * s))
    fac = max(1 - sy / se, 0.0) if se > 0 else 0.0
    return (1 + nu) / E * sig - nu / E * sk * np.eye(3) + 1.5 * (1 / T - 1 / E) * fac * s


def leading_compliance(T):
    """C[i,j,k,l] with eps_ij = C_ijkl sig_kl, from lim eps(L sig)/L; basis stresses symmetric."""
    L = mp.mpf(10) ** 25
    C = np.zeros((3, 3, 3, 3))
    for k in range(3):
        for l in range(k, 3):
            B = mp.zeros(3, 3); B[k, l] = 1; B[l, k] = 1
            e = full_law_mp(L * B, T) / L
            for i in range(3):
                for j in range(3):
                    v = float(e[i, j])
                    if k == l:
                        C[i, j, k, l] = v
                    else:
                        C[i, j, k, l] = C[i, j, l, k] = v / 2
    return C


def westergaard(x, y, K):
    z = complex(x, y)
    Z = K / np.sqrt(2 * np.pi * z)          # principal branch: cut along x < 0
    Zp = -0.5 * Z / z
    sxx = Z.real - y * Zp.imag
    syy = Z.real + y * Zp.imag
    sxy = -y * Zp.real
    return np.array([[sxx, sxy, 0.0], [sxy, syy, 0.0], [0.0, 0.0, 0.0]])


class Field:
    def __init__(self, T, K=1.0):
        self.T, self.K, self.C = T, K, leading_compliance(T)

    def sig(self, x, y):
        return westergaard(x, y, self.K)

    def eps(self, x, y):
        return np.einsum("ijkl,kl->ij", self.C, self.sig(x, y))

    def deps(self, x, y, d):
        return [(self.eps(x + d, y) - self.eps(x - d, y)) / (2 * d),
                (self.eps(x, y + d) - self.eps(x, y - d)) / (2 * d)]

    def J_half(self, h, upper=True):
        """J over the upper (or lower) half of the square |x|,|y| <= h, from the ligament point
        (h, 0) to the face; omega integrated along the path (omega = 0 on y = 0, x > 0)."""
        sgn = 1.0 if upper else -1.0
        legs = [((h, 0.0), (h, sgn * h), (1.0, 0.0)),
                ((h, sgn * h), (-h, sgn * h), (0.0, sgn)),
                ((-h, sgn * h), (-h, sgn * 1e-13), (-1.0, 0.0))]
        omega, Jtot, d = 0.0, 0.0, 1e-6 * h
        for (x0, y0), (x1, y1), N in legs:
            Ls = np.hypot(x1 - x0, y1 - y0); tx, ty = (x1 - x0) / Ls, (y1 - y0) / Ls

            def rhs(s, yv):
                x, y = x0 + tx * s, y0 + ty * s
                (ex, ey) = self.deps(x, y, d)
                om_x = ex[0, 1] - ey[0, 0]           # d omega/dx = eps_xy,x - eps_xx,y
                om_y = ex[1, 1] - ey[0, 1]           # d omega/dy = eps_yy,x - eps_xy,y
                sg, ep = self.sig(x, y), self.eps(x, y)
                W = 0.5 * np.sum(sg * ep)
                t = sg[:2, :2] @ np.array(N)
                uxx = ep[0, 0]; uyx = ep[0, 1] + yv[0]
                return [om_x * tx + om_y * ty, W * N[0] - (t[0] * uxx + t[1] * uyx)]
            sol = solve_ivp(rhs, (0, Ls), [omega, 0.0], method="DOP853", rtol=1e-11, atol=1e-16)
            omega = sol.y[0, -1]; Jtot += sol.y[1, -1]
        return Jtot

    def J(self, h):
        return self.J_half(h, True) + self.J_half(h, False)


def coeffs(T, K):
    f = Field(T, K)
    r = 1.0
    out = {}
    out["C_zz0"] = f.eps(r, 0.0)[2, 2] * np.sqrt(r)
    th = np.pi / 3
    R = np.array([[np.cos(th), np.sin(th), 0], [-np.sin(th), np.cos(th), 0], [0, 0, 1]])   # rows e_r, e_th
    ep = R @ f.eps(r * np.cos(th), r * np.sin(th)) @ R.T
    out["C_rr"], out["C_tt"], out["C_g"] = ep[0, 0], ep[1, 1], 2 * ep[0, 1]
    yv = 0.37
    uy = quad(lambda t: 2 * t * f.eps(0.0, t * t)[1, 1], 0.0, np.sqrt(yv), epsabs=0, epsrel=1e-13)[0]
    out["C_y"] = uy / np.sqrt(yv)
    # screens on theta = 0
    se = lambda rr: (lambda s: np.sqrt(1.5 * np.sum((s - np.trace(s) / 3 * np.eye(3)) ** 2)))(f.sig(rr, 0.0))
    out["r_max_stress"] = brentq(lambda rr: se(rr) - 10 * sy, 1e-8, 10.0, xtol=1e-16, rtol=1e-14)
    pmax = lambda rr: np.max(np.abs(np.linalg.eigvalsh(f.eps(rr, 0.0))))
    out["r_min_strain"] = brentq(lambda rr: pmax(rr) - 0.05, 1e-8, 10.0, xtol=1e-16, rtol=1e-14)
    out["r_max"] = min(out["r_max_stress"], 0.200)
    out["r_min"] = out["r_min_strain"]
    return out


if __name__ == "__main__":
    res = {}
    # 1. leading compliance: isotropic? linear? and its uniaxial slope
    T0 = 5000.0
    C0 = leading_compliance(T0)
    s1 = np.diag([1.0, 0, 0])
    e1 = np.einsum("ijkl,kl->ij", C0, s1)
    res["uniaxial_slope_check_T5000"] = {"1/eps11": 1 / e1[0, 0], "-eps22/eps11": -e1[1, 1] / e1[0, 0]}
    # 2. J per K^2 on two square contours, at T0 (numeric), and the implied relation
    f0 = Field(T0, 1.0)
    j1, j2 = f0.J(1.0), f0.J(0.3)
    res["J_per_K2_T5000"] = {"square_h1": j1, "square_h0.3": j2, "times_T": j1 * T0}
    # 3. inverse: C_zz(T) with K from the numeric J; scan for monotonicity, then brentq
    def czz_of_T(T):
        f = Field(T, 1.0)
        jk = f.J(0.5)
        K = np.sqrt(Jgiven / jk)
        return f.eps(1.0, 0.0)[2, 2] * K, K
    grid = np.geomspace(10.0, E * 0.999, 25)
    vals = [czz_of_T(T)[0] for T in grid]
    res["Czz_scan"] = {"T": list(grid), "Czz": vals,
                       "monotone_increasing": bool(np.all(np.diff(vals) > 0))}
    T = brentq(lambda T: czz_of_T(T)[0] - Czz_given, 100.0, E * 0.999, xtol=1e-12, rtol=1e-15)
    K = czz_of_T(T)[1]
    out = coeffs(T, K)
    out["T"] = T; out["K"] = K
    res["values"] = out
    key = {"T": 4816.472187302739, "C_rr": 0.003942048237733179, "C_tt": 0.0006331781350710297,
           "C_g": 0.005731131133456488, "C_y": 0.009139018941480233, "r_min": 0.009999999999999998,
           "r_max": 0.06132522854990914}   # reference_output.json route A
    res["rel_diff_vs_key"] = {k: out[k] / key[k] - 1 for k in key}
    # 4. closure: J of the final field on a third contour
    res["J_final_square_h0.07"] = Field(T, K).J(0.07)
    json.dump(res, open(os.path.join(HERE, "recompute_V2_results.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "Czz_scan"}, indent=1))
    print("Czz monotone on (10, E):", res["Czz_scan"]["monotone_increasing"],
          "range", vals[0], vals[-1])
