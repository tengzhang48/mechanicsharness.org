"""Sharp V-notch (opening 120 deg) in plane-strain power-law material.

Generalizes the HRR crack analysis (knowledge/hrr.py, Hutchinson 1968) to a notch:
material occupies -alpha < theta < alpha with alpha = 2*pi/3, traction-free notch
flanks at theta = +-alpha, symmetric (mode-I-type) remote tension.

Airy  phi = K r^s f(theta):
    sr  = r^(s-2) (s f + f'')
    st  = r^(s-2) s(s-1) f
    srt = r^(s-2) (1-s) f'
Plane strain, plastic (incompressible-deviatoric) flow dominates as r->0:
    se^2 = (3/4)(sr-st)^2 + 3 srt^2,   eps_ij ~ se^(n-1) s_ij (nonlinear)
ODE (Hutchinson eq. 31 with a = n(s-2), b = a+1):
    [d2/dth2 - a(a+2)] [ se^(n-1) (s(2-s) f + f'') ]
        + 4 (s-1) b ( se^(n-1) f' )' = 0
BCs: symmetry at 0: f'(0) = f'''(0) = 0; traction-free at alpha: f(alpha) = f'(alpha) = 0.
Unknowns: s (eigenvalue; p = 2 - s) and f''(0), with f(0) = 1 normalization.
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve

sys.path.insert(0, "knowledge")
import hrr

ALPHA = 2.0 * np.pi / 3.0
N = 4.0


def integrate(n, s, p0, alpha, rtol=1e-12, atol=1e-14, max_step=0.01):
    if abs(n - 1.0) < 1e-14:
        # closed linear limit: rho=1, g=0, R0=0 (hrr rhs gives 0*inf=nan at se=0)
        def rhs(x0, x1, x2, x3):
            C = s * (2.0 - s) * x0 + x2
            rest = (s * (2.0 - s) * x2 - (s - 2.0) * s * C
                    + 4.0 * (s - 1.0) ** 2 * x2)
            return -rest
    else:
        rhs = hrr._make_rhs(True, n, s)

    def f(t, y):
        with np.errstate(all="ignore"):
            v = rhs(y[0], y[1], y[2], y[3])
        return [y[1], y[2], y[3], v]
    return solve_ivp(f, [0.0, alpha], [1.0, 0.0, p0, 0.0],
                     rtol=rtol, atol=atol, dense_output=True, max_step=max_step)


def residual(x, n, alpha, **kw):
    s, p0 = x
    try:
        with np.errstate(all="ignore"):
            sol = integrate(n, s, p0, alpha, **kw)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
            return [1e6, 1e6]
        return [sol.y[0, -1], sol.y[1, -1]]
    except Exception:
        return [1e6, 1e6]


def elastic_symmetric_eigenvalue(alpha):
    """n->1 limit: characteristic equation
    (lam+1) tan((lam+1) a) = (lam-1) tan((lam-1) a); smallest lam in (1/2, 1).
    For alpha=pi the root is exactly lam=1/2 (endpoint of the interval)."""
    def g(lam):
        return (lam + 1.0) * np.tan((lam + 1.0) * alpha) \
             - (lam - 1.0) * np.tan((lam - 1.0) * alpha)
    if abs(alpha - np.pi) < 1e-14:
        return 0.5
    grid = np.linspace(0.5000001, 0.99999, 20001)
    vals = g(grid)
    idx = np.where(vals[:-1] * vals[1:] < 0)[0]
    if len(idx) == 0:
        raise RuntimeError("no sign change")
    i = idx[0]
    return brentq(g, grid[i], grid[i + 1], xtol=1e-15, rtol=8.9e-16)


def elastic_start(alpha):
    if abs(alpha - np.pi) < 1e-14:
        # crack: f = cos(th/2) + (1/3) cos(3th/2)  =>  f''(0)/f(0) = -3/4
        return 1.5, -0.75
    lam = elastic_symmetric_eigenvalue(alpha)
    C = -np.cos((lam + 1.0) * alpha) / np.cos((lam - 1.0) * alpha)
    p0 = (-(lam + 1.0) ** 2 - C * (lam - 1.0) ** 2) / (1.0 + C)
    return lam + 1.0, p0, C


def solve_notch(n, alpha, x0=None, **kw):
    if x0 is None:
        s0, p0, _ = elastic_start(alpha)
        x0 = [s0, p0]
    x, info, ier, msg = fsolve(lambda x: residual(x, n, alpha, **kw), x0,
                               full_output=True, xtol=1e-13)
    return x, ier, msg


def stresses(n, s, p0, alpha, n_th=4001, **kw):
    sol = integrate(n, s, p0, alpha, **kw)
    th = np.linspace(0.0, alpha, n_th)
    Y = sol.sol(th)
    f, fp, fpp = Y[0], Y[1], Y[2]
    sr = s * f + fpp
    st = s * (s - 1.0) * f
    srt = (1.0 - s) * fp
    se = np.sqrt(0.75 * (sr - st) ** 2 + 3.0 * srt ** 2)
    return dict(theta=th, f=f, fp=fp, fpp=fpp, sr=sr, st=st, srt=srt, se=se,
                s=s, p0=p0, sol=sol)


def main():
    alpha = ALPHA

    # --- start: n=1 limit matches analytic characteristic root
    lam = elastic_symmetric_eigenvalue(alpha)
    print(f"elastic(incompressible-limit) symmetric eigenvalue: lam={lam:.12f}  p=1-lam={1-lam:.12f}")

    # --- continuation in n from 1 to 4
    n_grid = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    x = None
    for n in n_grid:
        if x is None:
            s0, p0, _ = elastic_start(alpha)
            x = [s0, p0]
        x, ier, msg = solve_notch(n, alpha, x0=x)
        s, p0 = x
        print(f"n={n:5.2f}  s={s:.12f}  p=2-s={2-s:.12f}  f''(0)={p0:.12f}  ier={ier}")

    n = N
    s, p0 = x
    p = 2.0 - s
    R = (s * 1.0 + p0) / (s * (s - 1.0))
    print()
    print(f"alpha = 2pi/3, n = {n}")
    print(f"p     = {p:.10f}")
    print(f"R     = {R:.10f}")
    print(f"check (n+1)p < 2 (finite energy): (n+1)p = {(n+1)*p:.6f}")

    # --- verification 1: crack limit alpha = pi must give s=(2n+1)/(n+1)
    lam_pi = elastic_symmetric_eigenvalue(np.pi)
    print(f"\n[check] elastic limit at alpha=pi: lam={lam_pi:.12f} (expect 0.5)")
    xc = list(elastic_start(np.pi)[:2])
    for nc in n_grid:
        xc, ier, _ = solve_notch(nc, np.pi, x0=xc)
    sc, pc0c = xc
    Rc = (sc + pc0c) / (sc * (sc - 1.0))
    print(f"[check] crack alpha=pi, n=4: s={sc:.12f} expect {hrr.singular_exponent(4.0):.12f}"
          f"  (p={2-sc:.12f} vs 1/(n+1)={1/5:.12f}); R_crack={Rc:.6f}")
    # n=3 crack angular-stress check against knowledge unit (Hutchinson Fig. 5a)
    x3, _ = solve_notch(3.0, np.pi, x0=[hrr.singular_exponent(3.0), -1.3])
    s3, p03 = x3
    print(f"[check] crack n=3: st(0)={s3*(s3-1):.4f} (expect 1.9407), "
          f"sr(0)={s3+p03:.4f} (expect ~1.65), ratio={((s3+p03)/(s3*(s3-1))):.4f}")

    # --- verification 2: residuals of the two notch BCs at convergence
    res = residual([s, p0], n, alpha)
    print(f"[check] notch BC residuals f(alpha), f'(alpha): {res[0]:.3e}, {res[1]:.3e}")

    # --- verification 2b: se > 0 on (0, alpha), opening stress positive
    d = stresses(n, s, p0, alpha)
    print(f"[check] min se on (0,alpha) = {d['se'].min():.4f} (>0); "
          f"st(0)/norm = {d['st'][0]:.4f} > 0; traction at flank: "
          f"st(alpha)={d['st'][-1]:.2e}, srt(alpha)={d['srt'][-1]:.2e}")

    # --- verification 2c: multistart search for other (more singular) branches
    found = []
    for s_guess in [1.65, 1.75, 1.84, 1.90, 1.95]:
        for p0_guess in [-1.0, -1.6, -2.2]:
            try:
                xx, ier, _ = fsolve(lambda z: residual(z, n, alpha, max_step=0.02),
                                    [s_guess, p0_guess], full_output=True,
                                    xtol=1e-12, maxfev=60)
            except Exception:
                continue
            if ier == 1 and max(abs(r) for r in residual(xx, n, alpha)) < 1e-7 \
                    and 1.0 < xx[0] < 2.0:
                if not any(abs(xx[0] - f0) < 1e-5 for f0, _ in found):
                    found.append((xx[0], xx[1]))
    found.sort()
    print(f"[check] eigenvalues found by multistart (s, R): "
          f"{[(round(sx, 8), round((sx + px) / (sx * (sx - 1)), 4)) for sx, px in found]}")

    # --- verification 3: ODE refinements (tolerance and max_step)
    for ms, rt, at in [(0.01, 1e-12, 1e-14), (0.002, 1e-13, 1e-15), (0.05, 1e-11, 1e-13)]:
        xx, ier, _ = solve_notch(n, alpha, x0=[s, p0], max_step=ms, rtol=rt, atol=at)
        pchk = 2 - xx[0]
        Rchk = (xx[0] + xx[1]) / (xx[0] * (xx[0] - 1.0))
        print(f"[check] max_step={ms}, rtol={rt:.0e}: p={pchk:.10f} R={Rchk:.10f}"
              f"  (dp={pchk-p:+.2e}, dR={Rchk-R:+.2e})")


if __name__ == "__main__":
    main()
