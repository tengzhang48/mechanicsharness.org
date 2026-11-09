"""Verification v3 (final): independent reconstruction + J, H, C, eigenvalue p.

Independent elements relative to knowledge/hrr.py:
- u_theta obtained by quadrature of u_theta' (eps_theta relation) from theta=0
  (odd symmetry), NOT the code's shear-strain formula.
- u_r' obtained analytically from dense-output states (no numeric diff).
- Cartesian gradients d/dx of u done two ways: my own closed-form polar->Cart
  formulas and central finite differences of u(x,y).
- Eigenvalue s (hence p) re-derived by two-parameter shooting (s, f''(0)).
"""
import sys, numpy as np
from scipy.optimize import brentq
sys.path.insert(0, "knowledge")
import hrr

alpha, sig0, Emod = 0.02, 400.0e6, 200.0e9
eps0 = sig0 / Emod

def make_case(n):
    p = 1.0 / (n + 1.0)
    s = hrr.singular_exponent(n)

    # --- two-parameter shooting: find s such that f(pi)=0 AND f'(pi)=0 -------
    def shoot_fpp0(s_try):
        try:
            return brentq(lambda q: hrr._integrate(True, n, q, s=s_try).y[0, -1],
                          -2.0, -0.2, xtol=1e-13, rtol=1e-14)
        except Exception:
            return np.nan
    def g(s_try):
        q = shoot_fpp0(s_try)
        if not np.isfinite(q):
            return np.nan
        return hrr._integrate(True, n, q, s=s_try).y[1, -1]
    s_grid = np.linspace(1.55, 1.95, 41)
    gv = np.array([g(x) for x in s_grid])
    brackets = [(s_grid[i], s_grid[i+1]) for i in range(len(s_grid)-1)
                if np.isfinite(gv[i]) and np.isfinite(gv[i+1]) and gv[i]*gv[i+1] < 0]
    s_star = brentq(g, *brackets[0], xtol=1e-12) if brackets else np.nan
    p_star = 2.0 - s_star     # since s = 2 - 1/(n+1)
    print(f"n={n}: eigenvalue shooting -> s* = {s_star:.10f} "
          f"(exact {(2*n+1)/(n+1):.10f}), p* = 2-s* = {p_star:.10f} "
          f"vs 1/(n+1) = {1/(n+1):.10f}, diff = {p_star - 1/(n+1):+.2e}")

    # --- angular solution (dense output) ------------------------------------
    pshoot, sol_ivp = hrr._shoot(True, n)
    def fields(theta):
        Y = sol_ivp.sol(np.abs(theta))
        sg = np.where(theta < 0, -1.0, 1.0)
        Y = Y * np.vstack([np.ones_like(sg), sg, np.ones_like(sg), sg])
        f, fp, fpp, fppp = Y
        sr = s * f + fpp
        st = s * (s - 1.0) * f
        srt = (1.0 - s) * fp
        se = np.sqrt(0.75 * (sr - st) ** 2 + 3.0 * srt ** 2)
        return f, fp, fpp, fppp, sr, st, srt, se
    thf = np.linspace(0.0, np.pi, 400001)
    c = 1.0 / fields(thf)[7].max()
    def F(theta):
        return [q * c for q in fields(theta)]

    # --- displacement angular functions (independent construction) ----------
    def ut_all(theta):
        f, fp, fpp, fppp, sr, st, srt, se = F(theta)
        rho = se ** (n - 1.0)
        ur = rho * 0.75 * (sr - st) / p
        utp = rho * 0.75 * (st - sr) - ur
        return ur, utp
    def build_ut(theta):
        """odd u_theta(theta) with u_theta(0)=0 via cumulative trapezoid"""
        pos = np.argsort(np.abs(theta))
        a = np.abs(theta)
        order = np.argsort(a)
        a_s = a[order]
        ur_s, utp_s = ut_all(np.clip(theta[order], 0, None) * 0 + a_s)
        ut_s = np.concatenate([[0.0], np.cumsum((utp_s[1:] + utp_s[:-1]) / 2 * np.diff(a_s))])
        ut = ut_s[np.argsort(order)] * np.sign(theta)
        return ut
    def urp_analytic(theta):
        f, fp, fpp, fppp, sr, st, srt, se = F(theta)
        d = sr - st
        dp = s * (2.0 - s) * fp + fppp
        se2 = 0.75 * d ** 2 + 3.0 * srt ** 2
        se2p = 1.5 * d * dp + 6.0 * srt * (1.0 - s) * fp**0  # placeholder fix below
        se2p = 1.5 * d * dp + 6.0 * srt * ((1.0 - s) * fpp)   # d(srt)/dtheta = (1-s) f''
        rho = se ** (n - 1.0)
        rhop = 0.5 * (n - 1.0) * rho * se2p / se2
        return (rhop * 0.75 * d + rho * 0.75 * dp) / p
    return p, s, F, build_ut, ut_all, urp_analytic

def run(n, label):
    p, s, F, build_ut, ut_all, urp_analytic = make_case(n)
    Jval, R = 5.0e-5, 1.0e-4
    sol_code = hrr.solve_angular(n, plane_strain=True, n_theta=32001)
    I_code = hrr.integral_I(n, plane_strain=True, sol=sol_code)
    En_code = hrr.energy_coefficient_En(n, plane_strain=True, sol=sol_code, I=I_code)
    Amp = sig0 * (Jval / (alpha * eps0 * sig0 * I_code)) ** p

    # compatibility residual of the shear-strain relation (uses my u_theta)
    thc = np.linspace(-np.pi, np.pi, 64001)
    ur_c, utp_c = ut_all(thc)
    ut_c = build_ut(thc)
    f, fp, fpp, fppp, sr, st, srt, se = F(thc)
    rho = se ** (n - 1.0)
    resid = ut_c * (1.0 - p) - (urp_analytic(thc) - 3.0 * rho * srt)
    print(f"[{label}] shear-compat residual: max={np.abs(resid).max():.3e} "
          f"(scale |ut|max={np.abs(ut_c).max():.2f})")

    scale = alpha * eps0 * (Amp / sig0) ** n
    def u_cart(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        ur_t, utp_t = ut_all(theta)
        ut_t = build_ut(theta)
        ur, ut = scale * r**p * ur_t, scale * r**p * ut_t
        return ur * np.cos(theta) - ut * np.sin(theta), ur * np.sin(theta) + ut * np.cos(theta)

    Nth = 400001
    thg = np.linspace(-np.pi, np.pi, Nth)
    ur_t, utp_t = ut_all(thg)
    ut_t = build_ut(thg)
    urp_t = urp_analytic(thg)
    f, fp, fpp, fppp, sr, st, srt, se = F(thg)
    ct, stt_ = np.cos(thg), np.sin(thg)

    # analytic Cartesian gradients (my own derivation)
    ruxx = p * (ur_t * ct - ut_t * stt_) * ct - stt_ * (urp_t * ct - ur_t * stt_ - utp_t * stt_ - ut_t * ct)
    ruyx = p * (ur_t * stt_ + ut_t * ct) * ct - stt_ * (urp_t * stt_ + ur_t * ct + utp_t * ct - ut_t * stt_)
    dux_dx = scale * R**(p - 1.0) * ruxx
    duy_dx = scale * R**(p - 1.0) * ruyx

    # FD Cartesian gradients
    xg, yg = R * ct, R * stt_
    dx = 1e-6 * R
    dux_dx_fd = (u_cart(xg + dx, yg)[0] - u_cart(xg - dx, yg)[0]) / (2 * dx)
    duy_dx_fd = (u_cart(xg + dx, yg)[1] - u_cart(xg - dx, yg)[1]) / (2 * dx)
    gdiff = max(np.abs(dux_dx - dux_dx_fd).max(), np.abs(duy_dx - duy_dx_fd).max())
    print(f"[{label}] analytic-vs-FD gradient max diff = {gdiff:.3e} "
          f"(grad scale {np.abs(duy_dx).max():.3e})")

    sig_rr = Amp * R**(-p) * sr
    sig_rt = Amp * R**(-p) * srt
    W = (n / (n + 1.0)) * alpha * eps0 * sig0 * (Amp * R**(-p) * se / sig0) ** (n + 1.0)
    tx = sig_rr * ct - sig_rt * stt_
    ty = sig_rr * stt_ + sig_rt * ct
    for tag, gx, gy in [("analytic", dux_dx, duy_dx), ("FD", dux_dx_fd, duy_dx_fd)]:
        Jd = np.trapezoid((W * ct - (tx * gx + ty * gy)) * R, thg)
        Id = Jd / (alpha * eps0 * sig0 * (Amp / sig0) ** (n + 1.0))
        print(f"[{label}] J({tag})/J = {Jd/Jval:.10f}   I_direct = {Id:.8f} "
              f"vs I_code {I_code:.8f} (rel {Id/I_code-1:+.2e})")

    H = Amp * R**(-p) * F(np.array([0.0]))[5][0] / (sig0 * (Jval / (alpha * eps0 * sig0 * R)) ** p)
    shape_full = np.trapezoid(F(thg)[7] ** (n + 1.0), thg)
    C = (n / (n + 1.0)) * shape_full / I_code
    print(f"[{label}] I={I_code:.8f} En={En_code:.6f} H={H:.8f} C={C:.8f} "
          f"En/pi={En_code/np.pi:.8f}")
    return I_code, H, C

run(6.5, "n=6.5")
run(3.0, "n=3  ")
