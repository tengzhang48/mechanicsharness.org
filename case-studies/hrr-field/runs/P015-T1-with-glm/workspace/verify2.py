"""Verification v2: dense-output angular functions, u_theta by own quadrature.

Independence from knowledge-code paths:
- f, f', f'', f''' evaluated from the solve_ivp dense output (no interpolation).
- displacements built from scratch: u_r from eps_rr, u_theta by quadrature of
  u_theta' = (3/2) rho s_th - u_r  (eps_theta relation), NOT the code's
  shear-strain formula (used only as a cross-check).
- J evaluated on the full circle with Cartesian displacement gradients by
  central differences of u(x,y).
"""
import sys, numpy as np
sys.path.insert(0, "knowledge")
import hrr

n = 6.5
p = 1.0 / (n + 1.0)
s = hrr.singular_exponent(n)

def build(n, plane_strain=True):
    pshoot, sol_ivp = hrr._shoot(plane_strain, n)
    M = hrr._quadratic_coeff(plane_strain, n, s)
    def raw(theta):
        z = np.abs(theta)
        Y = sol_ivp.sol(z)          # dense output: f, f', f'', f'''
        sgn = np.sign(theta)
        Y = Y * np.vstack([sgn**0, sgn**1, sgn**2, sgn**3])
        return Y
    def fields(theta):
        Y = raw(theta)
        f, fp, fpp, fppp = Y
        sr = s * f + fpp
        st = s * (s - 1.0) * f
        srt = (1.0 - s) * fp
        se = np.sqrt(0.75 * (sr - st) ** 2 + 3.0 * srt ** 2)
        return f, fp, fpp, fppp, sr, st, srt, se
    # normalisation: max se = 1 on a fine grid of [0, pi]
    thf = np.linspace(0.0, np.pi, 200001)
    sef = fields(thf)[7]
    c = 1.0 / sef.max()
    def fields_n(theta):
        f, fp, fpp, fppp, sr, st, srt, se = fields(theta)
        return f*c, fp*c, fpp*c, fppp*c, sr*c, st*c, srt*c, se*c
    return fields_n, pshoot

def utildes(theta, fields):
    """u_r~ and u_theta~ from eps_rr and eps_theta_theta (independent route)."""
    f, fp, fpp, fppp, sr, st, srt, se = fields(theta)
    rho = se ** (n - 1.0)
    b = p
    ur = rho * 0.75 * (sr - st) / b          # (3/2)(n+1) rho s_r-dev / ...
    utp = rho * 0.75 * (st - sr) - ur        # eps_theta relation
    ut = np.concatenate([[0.0], np.cumsum((utp[1:] + utp[:-1]) / 2.0 * np.diff(theta))])
    return ur, ut, utp

def run_case(n, plane_strain=True, label=""):
    global p, s
    p = 1.0 / (n + 1.0); s = hrr.singular_exponent(n)
    fields, pshoot = build(n, plane_strain)

    # material / loading
    alpha, sig0, E = 0.02, 400.0e6, 200.0e9
    eps0 = sig0 / E
    Jval = 5.0e-5
    R = 1.0e-4
    sol_code = hrr.solve_angular(n, plane_strain=plane_strain, n_theta=16001)
    I_code = hrr.integral_I(n, plane_strain=plane_strain, sol=sol_code)
    Amp = sig0 * (Jval / (alpha * eps0 * sig0 * I_code)) ** p

    # ---- compatibility cross-check of the two u_theta formulas -------------
    thc = np.linspace(-np.pi, np.pi, 32001)
    ur_c, ut_c, utp_c = utildes(thc, fields)
    f, fp, fpp, fppp, sr, st, srt, se = fields(thc)
    rho = se ** (n - 1.0)
    urp_c = np.gradient(ur_c, thc)
    ut_shear = (urp_c - 3.0 * rho * srt) / (1.0 - p)   # code's shear route
    dmax = np.max(np.abs(ut_c - ut_shear))
    print(f"[{label}] max|u_th(eps-route) - u_th(shear-route)| = {dmax:.3e}")

    # ---- J direct ----------------------------------------------------------
    def u_cart(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        ur_t, ut_t, _ = utildes(theta, fields)
        scale = alpha * eps0 * (Amp / sig0) ** n * r ** p
        ur, ut = scale * ur_t, scale * ut_t
        return ur * np.cos(theta) - ut * np.sin(theta), ur * np.sin(theta) + ut * np.cos(theta)

    Nth = 200001
    thg = np.linspace(-np.pi, np.pi, Nth)
    xg, yg = R * np.cos(thg), R * np.sin(thg)
    dx = 1e-6 * R
    uxp = (u_cart(xg + dx, yg)[0] - u_cart(xg - dx, yg)[0]) / (2 * dx)
    uyp = (u_cart(xg + dx, yg)[1] - u_cart(xg - dx, yg)[1]) / (2 * dx)
    f, fp, fpp, fppp, sr, st, srt, se = fields(thg)
    sig_rr = Amp * R ** (-p) * sr
    sig_rt = Amp * R ** (-p) * srt
    W = (n / (n + 1.0)) * alpha * eps0 * sig0 * (Amp * R ** (-p) * se / sig0) ** (n + 1.0)
    tx = sig_rr * np.cos(thg) - sig_rt * np.sin(thg)
    ty = sig_rr * np.sin(thg) + sig_rt * np.cos(thg)
    J_direct = np.trapezoid((W * np.cos(thg) - (tx * uxp + ty * uyp)) * R, thg)

    I_direct = J_direct / (alpha * eps0 * sig0 * (Amp / sig0) ** (n + 1.0))
    print(f"[{label}] I_code = {I_code:.8f}   I_direct = {I_direct:.8f}   "
          f"rel = {I_direct/I_code-1:+.3e}")
    print(f"[{label}] J_direct/J_given = {J_direct/Jval:.10f}")

    # ---- H and C from the reconstructed field -------------------------------
    f0, fp0, fpp0, fppp0, sr0, st0, srt0, se0 = fields(np.array([0.0]))
    H = (Amp * R ** (-p) * st0[0]) / (sig0 * (Jval / (alpha * eps0 * sig0 * R)) ** p)
    th2 = np.linspace(-np.pi, np.pi, 200001)
    f2, fp2, fpp2, fppp2, sr2, st2, srt2, se2 = fields(th2)
    shape_full = np.trapezoid(se2 ** (n + 1.0), th2)
    C = (n / (n + 1.0)) * shape_full / I_code
    En_code = hrr.energy_coefficient_En(n, plane_strain=plane_strain, sol=sol_code, I=I_code)
    print(f"[{label}] H = {H:.8f}   C = {C:.8f}   E_n/pi = {En_code/np.pi:.8f}")
    return H, C

print("=== n = 6.5, plane strain (target case) ===")
run_case(6.5, True, "n=6.5 ps")
print("=== n = 3, plane strain (known: I=5.507, En=1.227) ===")
run_case(3.0, True, "n=3 ps  ")
