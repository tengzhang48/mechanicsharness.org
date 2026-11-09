"""Independent end-to-end verification of J, H, C.

Build the FULL dimensional stress and displacement HRR field, then compute
  J = int_{-pi}^{pi} (W cos(theta) - t_i du_i/dx) r dtheta
by finite differences (Cartesian gradients), and also compute U(R) by direct
2-D integration of W.  Compare to the closed-form relations used in solve.py.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

n = 6.5
alpha = 0.02
sig0 = 400e6
E = 200e9
eps0 = sig0 / E
s = (2 * n + 1) / (n + 1)


def Mcoeff(s):
    M = np.zeros((3, 3))
    M[0, 0] = 0.75 * s * s * (2 - s) ** 2
    M[1, 1] = 3 * (1 - s) ** 2
    M[2, 2] = 0.75
    M[0, 2] = M[2, 0] = 0.75 * s * (2 - s)
    return M


def make_rhs(n, s):
    k = (n - 1.0) / 2.0
    a = n * (s - 2.0)
    b = a + 1.0
    M = Mcoeff(s)

    def rhs(x0, x1, x2, x3):
        z = np.array([x0, x1, x2]); v = np.array([x1, x2, x3])
        Q = z @ M @ z; dQ = 2.0 * (M @ z); dQv = dQ @ v; vTMv = v @ M @ v
        rho = Q ** k; g = k * Q ** (k - 1.0) * dQ
        g0, g1, g2 = g
        R0 = (k * Q ** (k - 1.0) * (2.0 * vTMv)
              + k * (k - 1.0) * Q ** (k - 2.0) * dQv ** 2 + g0 * x2 + g1 * x3)
        P = g @ v
        C = s * (2.0 - s) * x0 + x2; Cp = s * (2.0 - s) * x1 + x3
        coef = rho + g2 * C
        rest = (R0 * C + 2.0 * P * Cp + rho * s * (2.0 - s) * x2
                - a * (a + 2.0) * rho * C + 4.0 * (s - 1.0) * b * (P * x1 + rho * x2))
        return -rest / coef
    return rhs


def solve_field():
    rhs = make_rhs(n, s)
    ff = lambda t, y: [y[1], y[2], y[3], rhs(y[0], y[1], y[2], y[3])]
    # 2D Newton for (s, p)
    def fend(ss, pp):
        sol = solve_ivp(ff, [0.0, np.pi], [1.0, 0.0, pp, 0.0],
                        rtol=1e-11, atol=1e-13, dense_output=True, max_step=0.004)
        return sol.y[0, -1], sol.y[1, -1]
    x = np.array([s + 1e-4, -0.689256761930])
    for _ in range(40):
        F = np.array(fend(x[0], x[1]))
        Jm = np.zeros((2, 2))
        for j in range(2):
            dx = np.zeros(2); dx[j] = 1e-8
            Jm[:, j] = (np.array(fend(*(x + dx))) - F) / 1e-8
        try:
            step = np.linalg.solve(Jm, -F)
        except np.linalg.LinAlgError:
            break
        x = x + step
        if np.linalg.norm(step) < 1e-13:
            break
    sol = solve_ivp(ff, [0.0, np.pi], [1.0, 0.0, x[1], 0.0],
                    rtol=1e-12, atol=1e-14, dense_output=True, max_step=0.004)
    return x[0], x[1], sol


def angular(sol_ivp, th):
    th = np.atleast_1d(np.asarray(th, dtype=float))
    Y = sol_ivp.sol(th)
    f, fp, fpp, fppp = Y
    M = Mcoeff(s)
    se = np.sqrt(np.einsum('ij,jk,ik->i', Y[:3].T, M, Y[:3].T))
    c = 1.0 / se.max()
    Y = Y * c
    f, fp, fpp, fppp = Y
    sr = s * f + fpp
    st = s * (s - 1.0) * f
    srt = (1.0 - s) * fp
    se = np.sqrt(np.einsum('ij,jk,ik->i', Y[:3].T, M, Y[:3].T))
    # displacements (plane strain), same as integral_I
    m = n * (s - 2.0); b = n * (s - 2.0) + 1.0
    rho = se ** (n - 1.0)
    z = np.stack([f, fp, fpp], -1)
    dQ = 2.0 * (z @ M.T)
    dQv = np.einsum('ij,ij->i', dQ, np.stack([fp, fpp, fppp], -1))
    rhop = (n - 1.0) * rho * (dQv / (2.0 * se ** 2))
    c1 = 0.75 * (sr - st); c1p = 0.75 * (s * (2.0 - s) * fp + fppp)
    ur = rho * c1 / b; urp = (rhop * c1 + rho * c1p) / b
    ut = (3.0 * rho * srt - urp) / m
    utp = rho * (0.75 * (st - sr)) - ur
    return dict(f=f, fp=fp, fpp=fpp, fppp=fppp, sr=sr, st=st, srt=srt, se=se,
                ur=ur, urp=urp, ut=ut, utp=utp, c=c)


s_num, p_shoot, sol_ivp = solve_field()
print("s =", s_num, " f''(0) =", p_shoot)

# amplitude: choose Kdim so J = 1 J/m^2 conceptually; but for verification just
# pick Kdim and check the J relation.
# Angular I and shape:
th = np.linspace(0.0, np.pi, 80001)
A = angular(sol_ivp, th)
n_ = n
integrand = ((n / (n + 1.0)) * A['se'] ** (n + 1.0) * np.cos(th)
             - (np.sin(th) * (A['sr'] * (A['ut'] - A['urp']) - A['srt'] * (A['ur'] + A['utp']))
                + (n * (s - 2.0) + 1.0) * np.cos(th) * (A['sr'] * A['ur'] + A['srt'] * A['ut'])))
I = 2.0 * np.trapezoid(integrand, th)
shape2 = 2.0 * np.trapezoid(A['se'] ** (n + 1.0), th)
print("I =", I, " shape2 =", shape2)

# =====================================================================
# Compute the DIMENSIONLESS J integral D of the normalized angular field
# (max se_tilde = 1), using the same rotation/finite-difference machinery,
# and compare with the eq.(24) I from solve.py.
#   dimensionless constitutive: eps~ = alpha * se~^(n-1) * s~   (se~,s~ /sig0;
#   eps~ /eps0); dimensionless displacement u~ with eps~ = grad u~.
#   Field: sigma~ = r^(s-2) sigma_ang,  u~ = r^(n(s-2)+1) u_ang (analytic part).
#   J_dimless D = int [ W~ cos th - t~ . du~/dx~ ] r dth.
# =====================================================================
def D_from_angular(sol_ivp, R=1.0, nth=80001, h=1e-6):
    th = np.linspace(-np.pi, np.pi, nth)
    x = R * np.cos(th); y = R * np.sin(th)

    def ang_at(xx, yy):
        t = np.arctan2(yy, xx)
        tt = np.where(t < 0, t + 2 * np.pi, t)
        tt = np.where(tt > np.pi, 2 * np.pi - tt, tt)
        return tt

    def stress_dimless(xx, yy, rr):
        tt = ang_at(xx, yy)
        Aa = angular(sol_ivp, tt)
        sc = rr ** (s - 2.0)
        sr = sc * Aa['sr']; st = sc * Aa['st']; srt = sc * Aa['srt']
        ct = xx / rr; sy = yy / rr
        sxx = sr * ct ** 2 + st * sy ** 2 - 2 * srt * ct * sy
        syy = sr * sy ** 2 + st * ct ** 2 + 2 * srt * ct * sy
        sxy = (sr - st) * ct * sy + srt * (ct ** 2 - sy ** 2)
        return sxx, syy, sxy

    def disp_dimless(xx, yy):
        rr = np.hypot(xx, yy)
        tt = ang_at(xx, yy)
        Aa = angular(sol_ivp, tt)
        aexp = n * (s - 2.0) + 1.0
        sc = alpha * rr ** aexp          # alpha so W and traction work share it
        ur = sc * Aa['ur']; ut = sc * Aa['ut']
        ct = xx / rr; sy = yy / rr
        return ur * ct - ut * sy, ur * sy + ut * ct

    # gradients
    def ddx(fun):
        return (fun(x + h, y) - fun(x - h, y)) / (2 * h)

    def ddy(fun):
        return (fun(x, y + h) - fun(x, y - h)) / (2 * h)

    up = disp_dimless(x + h, y); um = disp_dimless(x - h, y)
    dux_dx = (up[0] - um[0]) / (2 * h)
    duy_dx = (up[1] - um[1]) / (2 * h)
    sxx, syy, sxy = stress_dimless(x, y, np.hypot(x, y))
    se = np.sqrt(0.75 * (sxx - syy) ** 2 + 3.0 * sxy ** 2)
    Wf = (n / (n + 1.0)) * alpha * se ** (n + 1.0)
    nx = np.cos(th); ny = np.sin(th)
    tx = sxx * nx + sxy * ny
    ty = sxy * nx + syy * ny
    integ = (Wf * np.cos(th) - (tx * dux_dx + ty * duy_dx)) * R
    return np.trapezoid(integ, th)

D = D_from_angular(sol_ivp)
print(f"Dimensionless J (numeric) D = {D:.8f}")
print(f"eq.(24) I (solve.py)        = {I:.8f}")
print(f"ratio I/D                   = {I/D:.8f}")


# =====================================================================
# FULL dimensional field and Cartesian finite-difference J integral.
# =====================================================================
p = 1.0 / (n + 1.0)

def Kdim_from_J(Jval):
    return (Jval * sig0 ** n / (alpha * eps0 * I)) ** (1.0 / (n + 1.0))

Kdim = 1.0e6   # arbitrary amplitude, stress*length^(2-s) units

def field_disp(x, y):
    """u_x, u_y of the dimensional HRR field at (x,y) [metres]."""
    r = np.hypot(x, y)
    t = np.arctan2(y, x)
    tt = np.where(t < 0, t + 2 * np.pi, t)          # map to [0,2pi)
    tt = np.where(tt > np.pi, 2 * np.pi - tt, tt)   # symmetry theta in [0,pi]
    Aa = angular(sol_ivp, tt)
    # amplitude: u_i = alpha eps0 (Kdim/sig0)^n r^(n(s-2)+1) utilde_i
    aexp = n * (s - 2.0) + 1.0
    scale = alpha * eps0 * (Kdim / sig0) ** n * r ** aexp
    ur = scale * Aa['ur']
    ut = scale * Aa['ut']
    ct = x / r; stt = y / r
    return ur * ct - ut * stt, ur * stt + ut * ct

def stress_cart(x, y):
    """sigma_xx, sigma_yy, sigma_xy of the dimensional HRR field."""
    r = np.hypot(x, y)
    t = np.arctan2(y, x)
    tt = np.where(t < 0, t + 2 * np.pi, t)
    tt = np.where(tt > np.pi, 2 * np.pi - tt, tt)
    Aa = angular(sol_ivp, tt)
    scale = Kdim * r ** (s - 2.0)
    sr = scale * Aa['sr']; st = scale * Aa['st']; srt = scale * Aa['srt']
    ct = x / r; stt = y / r; c2 = ct * ct; s2 = stt * stt; cs = ct * stt
    sxx = sr * c2 + st * s2 - 2 * srt * cs
    syy = sr * s2 + st * c2 + 2 * srt * cs
    sxy = (sr - st) * cs + srt * (c2 - s2)
    return sxx, syy, sxy

def eff_stress(x, y):
    sxx, syy, sxy = stress_cart(x, y)
    return np.sqrt(0.75 * (sxx - syy) ** 2 + 3.0 * sxy ** 2)

def W_energy(x, y):
    se = eff_stress(x, y)
    return (n / (n + 1.0)) * alpha * eps0 * se ** (n + 1.0) / sig0 ** n

def J_numeric_cartesian(R, nth=20001, h=1e-7):
    """Circular-contour J with finite-difference Cartesian gradients."""
    th = np.linspace(-np.pi, np.pi, nth)
    x = R * np.cos(th)
    y = R * np.sin(th)
    # displacement gradients du_i/dx via central differences in x
    ux = field_disp(x, y)
    ux_p = field_disp(x + h, y)
    ux_m = field_disp(x - h, y)
    dux_dx = (ux_p[0] - ux_m[0]) / (2 * h)
    duy_dx = (ux_p[1] - ux_m[1]) / (2 * h)
    # traction t_i = sigma_ij N_j with N = radial unit vector (cos, sin)
    sxx, syy, sxy = stress_cart(x, y)
    nx = np.cos(th); ny = np.sin(th)
    tx = sxx * nx + sxy * ny
    ty = sxy * nx + syy * ny
    W = W_energy(x, y)
    integ = (W * np.cos(th) - (tx * dux_dx + ty * duy_dx)) * R
    return np.trapezoid(integ, th)

for R in [1e-4, 1e-3, 1e-2]:
    Jn = J_numeric_cartesian(R)
    Jclosed = alpha * eps0 / sig0 ** n * Kdim ** (n + 1.0) * I
    print(f"R={R:.0e}: J_numeric={Jn:.6e}  J_closed={Jclosed:.6e}  ratio={Jn/Jclosed:.8f}")