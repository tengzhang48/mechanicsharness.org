"""HRR crack-tip field and notch variant (Hutchinson 1968, corpus P015): recompute the reference numbers.

Written on 25 October 2026 from the problem statements (probe tasks T1, T2; matched tests V1, V2) and the
paper's formulation (plane-strain angular equation, Eq. 31; J contour, Eqs. 23-24), without the builders'
reference.py files or the benchmark audit's recompute_key.py. It checks:
  1. T1: plane-strain crack, fully plastic J2 power law, n = 6.5: exponent with s left free, H, C, I_n, E_n;
     the n = 1 limit; Table 1 of the paper (plane strain, n = 3, 5, 9, 13); a compatibility residual.
  2. T2: 120-degree notch (material in |theta| < 2pi/3), n = 4: p and R; limits n = 1 (elastic wedge root),
     faces at +-pi (crack, p = 1/(n+1)); the paper's 90-degree notch at n = 3 (printed exponent 0.225).
  3. Curves for the page figure: p(n) for the crack and the 120- and 90-degree notches, and I_n(n), n = 1..13.
  4. V1: antiplane power-law crack, n = 3.5, J = 0.120 N/mm: the seven requested numbers, with J evaluated in
     Cartesian components on a circle (a mixed-variable shooting solve, not the builder's routes A or B).
  5. V2: bilinear (finite tangent stiffness) membrane: reduction of the leading compliance, inverse for T, the
     six other requested numbers.
  6. Run configuration: event types in the DeepSeek (Cline) session logs of the probe round and reasoning-token
     totals in the DeepSeek/GLM (OpenCode) logs of the matched test. Counts only; no log content is printed.
Run:  python -I hrr_check.py > check-outputs/hrr_check.txt
"""
import json
import os
import platform
import re
import time

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, fsolve

RT = dict(method="DOP853", rtol=1e-12, atol=1e-13)


# ---------------------------------------------------------------- in-plane (plane strain) power-law field
def invert_D(g, v, h, n):
    """Solve Q^(n-1) D = g for D, with Q^2 = 3/4 D^2 + 3 h^2 v^2 (monotone in D)."""
    if n == 1.0 or g == 0.0:
        return g
    c = 3.0 * h * h * v * v
    a = abs(g)
    d = a ** (1.0 / n) / 0.75 ** ((n - 1) / (2 * n))           # bound from Q >= sqrt(3/4)|D|
    if c > 0:
        d = min(d, a / c ** ((n - 1) / 2))                       # bound from Q >= sqrt(c)
    for _ in range(200):                                         # Newton from above (convex, increasing)
        Q2 = 0.75 * d * d + c
        psi = Q2 ** ((n - 1) / 2) * d - a
        dpsi = Q2 ** ((n - 3) / 2) * (Q2 + (n - 1) * 0.75 * d * d)
        step = psi / dpsi
        d -= step
        if abs(step) <= 1e-15 * max(d, 1e-300):
            break
    return np.copysign(d, g)


def rhs_inplane(th, y, s, n):
    """y = (f, v = f', g = Q^(n-1) D, G = g' + 4h(k+1) m), m = Q^(n-1) f'; compatibility gives G' = k(k+2) g."""
    f, v, g, G = y
    h, k = s - 1.0, n * (s - 2.0)
    D = invert_D(g, v, h, n)
    Q = np.sqrt(0.75 * D * D + 3 * h * h * v * v)
    m = Q ** (n - 1) * v
    return [v, D - s * (2 - s) * f, G - 4 * h * (k + 1) * m, k * (k + 2) * g]


def g_of(fpp0, s, n):
    D0 = fpp0 + s * (2 - s)
    return (0.75 * D0 * D0) ** ((n - 1) / 2) * D0


def shoot(beta, n, s0, q0):
    """Unknowns s and q = f''(0); f(0) = 1, f'(0) = 0, symmetric (G(0) = 0); f(beta) = f'(beta) = 0."""
    def res(x):
        s, q = x
        sol = solve_ivp(rhs_inplane, (0, beta), [1.0, 0.0, g_of(q, s, n), 0.0], args=(s, n), **RT)
        return [sol.y[0, -1], sol.y[1, -1]]
    x, info, ier, msg = fsolve(res, [s0, q0], full_output=True, xtol=1e-14)
    resid = max(abs(r) for r in res(x))
    if resid > 1e-10:
        raise RuntimeError(f"no convergence beta={beta} n={n}: {msg}")
    return x[0], x[1], resid


def crack_numbers(n, s, q, npts=4001):
    """J coefficient J0 (f(0) = 1), I_n (max Q = 1), H, C, E_n = pi C, and a compatibility residual."""
    th = np.linspace(0, np.pi, npts)
    sol = solve_ivp(rhs_inplane, (0, np.pi), [1.0, 0.0, g_of(q, s, n), 0.0], args=(s, n),
                    dense_output=True, **RT)
    f, v, g, G = sol.sol(th)
    h, k = s - 1.0, n * (s - 2.0)
    a = k + 1.0                                                  # displacement exponent
    D = np.array([invert_D(gi, vi, h, n) for gi, vi in zip(g, v)])
    Q = np.sqrt(0.75 * D * D + 3 * h * h * v * v)
    m = Q ** (n - 1) * v
    gp = G - 4 * h * (k + 1) * m
    er, erth = 0.75 * g, -1.5 * h * m                            # eps_rr, eps_rtheta (tensor) angular parts
    Gu, Gup = er / a, 0.75 * gp / a                              # u_r = r^a Gu
    L = (2 * erth - Gup) / (a - 1.0)                             # u_theta = r^a L, from the shear relation
    Lp = -er - Gu                                                # from eps_thth = -eps_rr
    srr, srt = D + s * (s - 1) * f, -h * v
    W = n / (n + 1) * Q ** (n + 1)
    uxr = a * Gu * np.cos(th) - (Gup - L) * np.sin(th)
    uxt = a * L * np.cos(th) - (Lp + Gu) * np.sin(th)
    integrand = W * np.cos(th) - srr * uxr - srt * uxt
    J0 = 2 * np.trapezoid(integrand, th)                         # integrand even in theta
    Qmax = Q.max()
    In = J0 / Qmax ** (n + 1)
    H = s * (s - 1) / J0 ** (1 / (n + 1))
    C = 2 * np.trapezoid(W, th) / J0
    compat = np.max(np.abs(np.gradient(L, th)[5:-5] - Lp[5:-5])) / np.max(np.abs(Lp))
    return dict(J0=J0, I=In, H=H, C=C, En=np.pi * C, Qmax=Qmax, compat_residual=compat,
                sigma_tt0_maxQ1=s * (s - 1) / Qmax)


def notch_R(s, q):
    return (s + q) / (s * (s - 1))


# ---------------------------------------------------------------- report
def main():
    out = []
    P = out.append
    P(f"# hrr_check.py run {time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())} UTC; python {platform.python_version()}, "
      f"numpy {np.__version__}, scipy {scipy.__version__}")
    data = {}

    # 1. T1 -----------------------------------------------------------------------------------------
    n = 6.5
    s, q, r = shoot(np.pi, n, 1.86, -0.7)
    c = crack_numbers(n, s, q)
    P("\n1. T1: plane-strain crack, n = 6.5, s left free (key: p 0.133333333333, H 1.893569725213, C 0.457305783974)")
    P(f"   s = {s:.12f} (exact (2n+1)/(n+1) = {(2*n+1)/(n+1):.12f}), p = 2 - s = {2-s:.12f}, f''(0) = {q:.10f}, face residual {r:.1e}")
    P(f"   J0 (f(0)=1) = {c['J0']:.10f}; I_n (max Q = 1) = {c['I']:.8f}; sigma_tt(0) with max Q = 1: {c['sigma_tt0_maxQ1']:.6f}")
    P(f"   H = {c['H']:.10f}; C = {c['C']:.10f}; E_n = pi C = {c['En']:.8f}; compatibility residual {c['compat_residual']:.1e}")
    data["T1"] = {"p": 2 - s, "H": c["H"], "C": c["C"], "I": c["I"], "En": c["En"]}
    s1, q1, _ = shoot(np.pi, 1.0, 1.5, -0.75)
    c1 = crack_numbers(1.0, s1, q1)
    P(f"   n = 1 limit: p = {2-s1:.10f} (exact 0.5), H = {c1['H']:.10f} (exact sqrt(2/(3 pi)) = {np.sqrt(2/(3*np.pi)):.10f}), "
      f"C = {c1['C']:.10f} (exact 0.25), I_1 = {c1['I']:.8f} (2 pi = {2*np.pi:.8f})")
    P("   Table 1 of the paper, plane strain (printed I, E_n on p. 21):")
    printed = {3: (5.51, 1.23), 5: (5.01, 1.37), 9: (4.60, 1.48), 13: (4.40, 1.52)}
    guess = {3: (1.75, -1.0), 5: (1.833, -0.8), 9: (1.9, -0.6), 13: (1.93, -0.5)}
    tab = {}
    for nn, (Ip, Ep) in printed.items():
        ss, qq, rr = shoot(np.pi, float(nn), *guess[nn])
        cc = crack_numbers(float(nn), ss, qq)
        tab[nn] = {"I": cc["I"], "En": cc["En"], "I_printed": Ip, "En_printed": Ep}
        P(f"   n = {nn:2d}: s = {ss:.10f} (exact {(2*nn+1)/(nn+1):.10f}); I = {cc['I']:.4f} vs {Ip}; E_n = {cc['En']:.4f} vs {Ep}; "
          f"differences {100*(cc['I']-Ip)/Ip:+.2f}% and {100*(cc['En']-Ep)/Ep:+.2f}%")
    data["table1"] = tab

    # 2. T2 -----------------------------------------------------------------------------------------
    P("\n2. T2: notch, material in |theta| < 2 pi/3, n = 4 (key: p 0.156175260286, R 0.373262313312)")
    b120 = 2 * np.pi / 3
    s2, q2, r2 = shoot(b120, 4.0, 1.84, -1.26)
    P(f"   s = {s2:.12f}, p = {2-s2:.12f}, f''(0) = {q2:.10f}, R = {notch_R(s2, q2):.12f}, face residual {r2:.1e}")
    P(f"   energy inside radius r grows as r^(2-(n+1)p) = r^{2-5*(2-s2):.6f} (crack: r^1)")
    data["T2"] = {"p": 2 - s2, "R": notch_R(s2, q2)}
    lam = brentq(lambda L: np.sin(2 * b120 * L) + L * np.sin(2 * b120), 0.55, 0.7)
    s21, q21, _ = shoot(b120, 1.0, 1.6, -1.17)
    A_ = 1.0
    # analytic symmetric wedge: Phi ~ r^(lam+1)[A cos((lam+1)th) + B cos((lam-1)th)], f(beta) = 0 fixes B/A
    B_ = -A_ * np.cos((lam + 1) * b120) / np.cos((lam - 1) * b120)
    f0 = A_ + B_
    fpp0 = -(A_ * (lam + 1) ** 2 + B_ * (lam - 1) ** 2)
    sl = lam + 1
    R_an = (sl * f0 + fpp0) / (sl * (sl - 1) * f0)
    P(f"   n = 1: p = {2-s21:.12f} vs elastic-wedge root 1 - lambda = {1-lam:.12f}; R = {notch_R(s21, q21):.12f} vs analytic {R_an:.12f}")
    s2c, q2c, _ = shoot(np.pi, 4.0, 1.8, -0.66)
    P(f"   faces at +-pi, n = 4: p = {2-s2c:.12f} (crack 1/(n+1) = 0.2), R = {notch_R(s2c, q2c):.6f}")
    s90, q90, _ = shoot(3 * np.pi / 4, 3.0, 1.775, -1.05)
    P(f"   paper's 90-degree notch, n = 3: p = {2-s90:.6f} (printed 0.225, p. 31), R = {notch_R(s90, q90):.6f}")

    # 3. curves -------------------------------------------------------------------------------------
    ns = [1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0]
    curves = {"n": ns, "crack_p": [], "crack_I": [], "notch120_p": [], "notch90_p": []}
    sc, qc = 1.5, -0.75
    s120, q120 = s21, q21
    s90c, q90c = shoot(3 * np.pi / 4, 1.0, 1.54, -1.1)[:2]
    for nn in ns:
        sc, qc, _ = shoot(np.pi, nn, sc, qc)
        curves["crack_p"].append(2 - sc)
        curves["crack_I"].append(crack_numbers(nn, sc, qc, npts=2001)["I"])
        s120, q120, _ = shoot(b120, nn, s120, q120)
        curves["notch120_p"].append(2 - s120)
        s90c, q90c, _ = shoot(3 * np.pi / 4, nn, s90c, q90c)
        curves["notch90_p"].append(2 - s90c)
    P("\n3. Curves for the figure (plane strain; crack p checked against 1/(n+1))")
    P("   n      crack p   1/(n+1)   I_n      notch120 p  notch90 p")
    for i, nn in enumerate(ns):
        P(f"   {nn:5.2f}  {curves['crack_p'][i]:.6f}  {1/(nn+1):.6f}  {curves['crack_I'][i]:.4f}   {curves['notch120_p'][i]:.6f}    {curves['notch90_p'][i]:.6f}")
    data["curves"] = curves

    # 4. V1 -----------------------------------------------------------------------------------------
    P("\n4. V1: antiplane power-law crack, n = 3.5, tau0 = 80 MPa, gamma0 = 0.002, J = 0.120 N/mm")
    n, tau0, gam0, Jv = 3.5, 80.0, 0.002, 0.120
    lam = 1 / (n + 1)
    pp = 1 + 1 / n
    kap = tau0 * gam0 ** (-1 / n)                                 # tau = kap |grad w|^(p-2) grad w

    def inv_fp(t, f):                                            # t = q^(p-2) f', q^2 = lam^2 f^2 + f'^2
        if t == 0.0:
            return 0.0
        if f == 0.0:                                             # q = |f'|: t = |f'|^(p-1) sign f'
            return np.copysign(abs(t) ** n, t)
        fun = lambda x: (lam * lam * f * f + x * x) ** ((pp - 2) / 2) * x - abs(t)   # increasing in x >= 0
        hi = 1.0
        while fun(hi) < 0:
            hi *= 2.0
        return np.copysign(brentq(fun, 0.0, hi, xtol=1e-16, rtol=1e-15), t)

    def rhs_ap(th, y):                                           # w = a r^lam f(theta), f odd; t ~ tau_theta
        f, t = y
        fp = inv_fp(t, f)
        qq = np.sqrt(lam * lam * f * f + fp * fp)
        return [fp, -lam * (1 - lam) * qq ** (pp - 2) * f]
    t0 = 1.0                                                     # f(0) = 0, f'(0) = 1 -> t(0) = 1
    sol = solve_ivp(rhs_ap, (0, np.pi), [0.0, t0], dense_output=True, method="DOP853", rtol=1e-12, atol=1e-14)
    th = np.linspace(0, np.pi, 8001)
    f, t = sol.sol(th)
    fp = np.array([inv_fp(ti, fi) for ti, fi in zip(t, f)])
    P(f"   face traction residual t(pi) = {t[-1]:.1e} (f'(0) = 1 normalization)")
    # J on the unit circle in Cartesian components, amplitude a = 1: gamma_r = lam f, gamma_th = f'
    gr, gt = lam * f, fp
    qq = np.sqrt(gr * gr + gt * gt)
    taur = kap * qq ** (pp - 2) * gr
    Wd = kap / pp * qq ** pp
    wx = gr * np.cos(th) - gt * np.sin(th)
    I_ap = 2 * np.trapezoid(Wd * np.cos(th) - taur * wx, th) / kap   # even integrand (f odd, f' even)
    a = (Jv / (kap * I_ap)) ** (1 / pp)
    P(f"   I = {I_ap:.12f} (closed form pi(n+1)/(2n) = {np.pi*(n+1)/(2*n):.12f}); a = {a:.12e} mm^(1-lam)")
    vals = {}
    fpi = sol.sol(np.pi)[0]
    for r_ in (0.02, 0.08):
        vals[f"slip_{r_}"] = 2 * a * r_ ** lam * fpi
    th60 = np.pi / 3
    f60, t60 = sol.sol(th60)
    fp60 = inv_fp(t60, f60)
    g = a * 0.02 ** (lam - 1) * np.array([lam * f60, fp60])
    tau = kap * np.hypot(*g) ** (pp - 2) * g
    vals["tau_rz"], vals["tau_tz"] = tau
    for ang, name in ((0.0, "r_0"), (np.pi / 2, "r_90"), (np.pi, "r_180")):
        fa, ta = sol.sol(ang)
        qa = np.hypot(lam * fa, inv_fp(ta, fa))
        vals[name] = (a * qa / 0.01) ** (1 / (1 - lam))
    key = {"slip_0.02": 0.001318126778, "slip_0.08": 0.001793693738, "tau_rz": 38.58659686,
           "tau_tz": 138.0106166, "r_0": 0.04689442305, "r_90": 0.02506612349, "r_180": 0.01339840659}
    for kname in key:
        P(f"   {kname:9s} = {vals[kname]:.10g}  key {key[kname]:.10g}  relative difference {abs(vals[kname]/key[kname]-1):.1e}")
    data["V1"] = vals

    # 5. V2 -----------------------------------------------------------------------------------------
    P("\n5. V2: membrane (sigma_zz = 0), complete law with finite tangent slope T; E = 70000, nu = 0.30, J = 0.800, C_zz(0) = -0.005")
    E, nu, Jm, Czz = 70000.0, 0.30, 0.800, -0.005

    def compliance(T, big=1e9):                                  # leading part of the stated tensor law
        def eps(sig):
            sk = np.trace(sig)
            dev = sig - sk / 3 * np.eye(3)
            return (1 + nu) / E * sig - nu / E * sk * np.eye(3) + 1.5 * (1 / T - 1 / E) * dev
        return eps
    nut = lambda T: 0.5 - (0.5 - nu) * T / E
    T_try = 5000.0
    sig = np.diag([3.0, -1.0, 0.5])
    e_full = compliance(T_try)(sig)
    e_iso = (1 + nut(T_try)) / T_try * sig - nut(T_try) / T_try * np.trace(sig) * np.eye(3)
    P(f"   leading compliance = isotropic (T, nu_t), nu_t = 1/2 - (1/2 - nu) T/E: max difference {np.abs(e_full-e_iso).max():.1e}")
    czz = lambda T: -2 * nut(T) * np.sqrt(Jm / (2 * np.pi * T))
    T = brentq(lambda T: czz(T) - Czz, 1.0, E, xtol=1e-12)
    nt = nut(T)
    K = np.sqrt(Jm * T)
    c, sn = np.cos(np.pi / 6), np.sin(np.pi / 6)
    S = K / np.sqrt(2 * np.pi) * np.array([c * (1 + sn * sn), c ** 3, sn * c * c])   # (rr, tt, rt) at 60 deg
    Crr = (S[0] - nt * S[1]) / T
    Ctt = (S[1] - nt * S[0]) / T
    Cg = 2 * (1 + nt) * S[2] / T
    mu, kap_t = T / (2 * (1 + nt)), (3 - nt) / (1 + nt)
    Cy = K / (2 * mu * np.sqrt(2 * np.pi)) * (kap_t - np.cos(np.pi / 2)) * np.sin(np.pi / 4)
    pr = np.sqrt(Jm / (2 * np.pi * T)) * np.array([1 - nt, 1 - nt, -2 * nt])        # principal strain coeffs, theta = 0
    rmin = (np.abs(pr).max() / 0.05) ** 2
    rmax = min(Jm * T / (2 * np.pi * 100.0 ** 2), 0.200)
    v2 = {"T": T, "C_rr": Crr, "C_tt": Ctt, "C_g": Cg, "C_y": Cy, "r_min": rmin, "r_max": rmax}
    key2 = {"T": 4816.472187, "C_rr": 0.003942048238, "C_tt": 0.0006331781351, "C_g": 0.005731131133,
            "C_y": 0.009139018941, "r_min": 0.01, "r_max": 0.06132522855}
    P(f"   nu_t = {nt:.10f}; largest principal strain on theta = 0 is the out-of-plane one: {abs(pr[2]) > abs(pr[0])}")
    for kname in key2:
        P(f"   {kname:6s} = {v2[kname]:.10g}  key {key2[kname]:.10g}  relative difference {abs(v2[kname]/key2[kname]-1):.1e}")
    P(f"   holding nu = 0.30 instead of nu_t would give T = {brentq(lambda T: -2*nu*np.sqrt(Jm/(2*np.pi*T)) - Czz, 1, E):.1f} MPa")
    data["V2"] = v2

    # 6. run configuration --------------------------------------------------------------------------
    P("\n6. Run configuration (counts only)")
    plog = "<probe-runs>-logs"
    for name in ("P015-T1-with-dsaz", "P015-T1-without-dsaz", "P015-T2-with-dsaz", "P015-T2-without-dsaz"):
        counts = {}
        with open(f"{plog}/{name}.jsonl", errors="replace") as fh:
            for line in fh:
                try:
                    d = json.loads(line)
                except ValueError:
                    continue
                e = d.get("event") if d.get("type") == "agent_event" else None
                if isinstance(e, dict) and e.get("contentType"):
                    counts[e["contentType"]] = counts.get(e["contentType"], 0) + 1
        P(f"   {name} (Cline): content events by type {counts}")
    llog = "<matched-runs>-logs"
    for name in ("P015-V1-B0-r1-k3az",):
        counts = {}
        with open(f"{llog}/{name}.jsonl", errors="replace") as fh:
            for line in fh:
                try:
                    d = json.loads(line)
                except ValueError:
                    continue
                e = d.get("event") if d.get("type") == "agent_event" else None
                if isinstance(e, dict) and e.get("contentType"):
                    counts[e["contentType"]] = counts.get(e["contentType"], 0) + 1
        P(f"   {name} (Cline, Kimi K3, for comparison): content events by type {counts}")
    for solver in ("deepseek", "glm"):
        tot = []
        for probe in ("V1", "V2"):
            for cond in ("B0", "B1", "C"):
                for rep in ("r1", "r2"):
                    name = f"P015-{probe}-{cond}-{rep}-{solver}"
                    rt = 0
                    with open(f"{llog}/{name}.jsonl", errors="replace") as fh:
                        for line in fh:
                            try:
                                d = json.loads(line)
                            except ValueError:
                                continue
                            if d.get("type") == "step_finish":
                                rt += ((d.get("part") or {}).get("tokens") or {}).get("reasoning", 0) or 0
                    tot.append(rt)
        P(f"   matched test, {solver} (OpenCode), 12 V1/V2 runs: reasoning tokens min {min(tot)}, max {max(tot)}")

    P("\nDATA " + json.dumps(data, default=float))
    print("\n".join(out))


if __name__ == "__main__":
    main()
