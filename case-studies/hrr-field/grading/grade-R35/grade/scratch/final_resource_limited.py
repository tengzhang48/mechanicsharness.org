"""FINAL, self-contained determination of p, H, C using ONLY numerical
integrals of the dimensional HRR field (validated against the elastic case).

Physical field (plane strain, total deformation law):
    sigma_ij(r,th) = sigma0 * A * r^(s-2) * s_tilde_ij(th)
    u_i(r,th)      = alpha*eps0 * A^n * r^(n(s-2)+1) * u_tilde_i(th)
with s = (2n+1)/(n+1), s_tilde normalised so max(se_tilde)=1, and u_tilde the
angular displacements consistent with the constitutive law + strain-displacement
relations (as used in integral_I of hrr.py).

J  = int_{-pi}^{pi} (W cos th - t.du/dx) r dth        (dimensional, J/m^2)
U  = int_0^R int W r dth dr                            (J/m)
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
p = 1.0 / (n + 1.0)


def Mcoeff(s):
    M = np.zeros((3, 3))
    M[0, 0] = 0.75 * s * s * (2 - s) ** 2
    M[1, 1] = 3 * (1 - s) ** 2
    M[2, 2] = 0.75
    M[0, 2] = M[2, 0] = 0.75 * s * (2 - s)
    return M


def make_rhs(n, s):
    k = (n - 1.0) / 2.0
    a = n * (s - 2.0); b = a + 1.0
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


rhs = make_rhs(n, s)
ff = lambda t, y: [y[1], y[2], y[3], rhs(y[0], y[1], y[2], y[3])]

# --- eigen value: s is EXACT analytically.  Fix s and solve for f''(0)=p by a
#     high-accuracy 1-D root find on f'(pi)=0 (equivalently f(pi)=0, both hold
#     at the eigenvalue).  This avoids the noisy 2-D FD Jacobian.
def yend(pp, tol=1e-12):
    sol = solve_ivp(ff, [0.0, np.pi], [1.0, 0.0, pp, 0.0],
                    rtol=tol, atol=1e-15, dense_output=False, max_step=0.002)
    return sol.y[0, -1], sol.y[1, -1]

# bracket on f'(pi)=0 near the known physical value
plo, phi = -0.70, -0.68
for _ in range(60):
    pm = 0.5 * (plo + phi)
    if yend(plo)[1] * yend(pm)[1] <= 0:
        phi = pm
    else:
        plo = pm
    if phi - plo < 1e-14:
        break
p_shoot = 0.5 * (plo + phi)
s_num = s
sol = solve_ivp(ff, [0.0, np.pi], [1.0, 0.0, p_shoot, 0.0],
                rtol=1e-12, atol=1e-15, dense_output=True, max_step=0.002)
print(f"eigen: s={s_num:.12f} (formula {s:.12f}), f''(0)={p_shoot:.12f}, "
      f"f(pi)={yend(p_shoot)[0]:.3e}, f'(pi)={yend(p_shoot)[1]:.3e}")


def _angular_raw(th):
    th = np.asarray(th, dtype=float)
    shape = th.shape
    tf = np.atleast_1d(th.ravel())
    Y = sol.sol(tf)
    return Y, shape

# Fixed normalisation constant: max se_tilde over a canonical grid [0,pi].
M = Mcoeff(s)
_Yc = sol.sol(np.linspace(0.0, np.pi, 20001))
_sec = np.sqrt(np.einsum('ij,jk,ik->i', _Yc[:3].T, M, _Yc[:3].T))
CNORM = 1.0 / _sec.max()


def angular(th):
    th = np.asarray(th, dtype=float)
    shape = th.shape
    tf = np.atleast_1d(th.ravel())
    Y = sol.sol(tf)
    se = np.sqrt(np.einsum('ij,jk,ik->i', Y[:3].T, M, Y[:3].T))
    c = CNORM                      # fixed, input-independent normalisation
    Y = Y * c
    f, fp, fpp, fppp = Y
    sr = s * f + fpp
    st = s * (s - 1.0) * f
    srt = (1.0 - s) * fp
    se = np.sqrt(np.einsum('ij,jk,ik->i', Y[:3].T, M, Y[:3].T))
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
    out = [a.reshape(shape) for a in (sr, st, srt, se, ur, urp, ut, utp)]
    return out


aexp = n * (s - 2.0) + 1.0     # = 1/(n+1) = p  (displacement radial exponent)


def fold(t):
    """Return (|th|, sign) with sign = +1 for th>=0, -1 for th<0.  Mode-I
    symmetry: sr, st, se, u_r are EVEN in th; srth, u_th are ODD.  The
    angular ODE is integrated on [0,pi] only, so continue by symmetry."""
    t = np.asarray(t, dtype=float)
    return np.abs(t), np.where(t < 0, -1.0, 1.0)


def sigma_cart(x, y, A):
    """Dimensional Cartesian stresses, amplitude A (metres^(2-s))."""
    r = np.hypot(x, y); t = np.arctan2(y, x)
    ta, sg = fold(t)
    sr, st, srt, se, *u = angular(ta)
    srt = srt * sg                     # odd component
    sc = sig0 * A * r ** (s - 2.0)
    SR, ST, SRT = sr * sc, st * sc, srt * sc
    ct = x / r; sy = y / r
    sxx = SR * ct ** 2 + ST * sy ** 2 - 2 * SRT * ct * sy
    syy = SR * sy ** 2 + ST * ct ** 2 + 2 * SRT * ct * sy
    sxy = (SR - ST) * ct * sy + SRT * (ct ** 2 - sy ** 2)
    return sxx, syy, sxy, SR, ST, SRT


def u_cart(x, y, A):
    r = np.hypot(x, y); t = np.arctan2(y, x)
    ta, sg = fold(t)
    *_, ur, urp, ut, utp = angular(ta)
    ut = ut * sg                       # odd component
    sc = alpha * eps0 * A ** n * r ** aexp
    UR, UT = ur * sc, ut * sc
    ct = x / r; sy = y / r
    return UR * ct - UT * sy, UR * sy + UT * ct


def se_cart(x, y, A):
    sxx, syy, sxy, *_ = sigma_cart(x, y, A)
    return np.sqrt(0.75 * (sxx - syy) ** 2 + 3.0 * sxy ** 2)


def W_cart(x, y, A):
    se = se_cart(x, y, A)
    return (n / (n + 1.0)) * alpha * eps0 * se ** (n + 1.0) / sig0 ** n


def J_numeric(A, R, nth=160001, h=1e-9):
    th = np.linspace(-np.pi, np.pi, nth)
    x = R * np.cos(th); y = R * np.sin(th)
    up = u_cart(x + h, y, A); um = u_cart(x - h, y, A)
    dux_dx = (up[0] - um[0]) / (2 * h)
    duy_dx = (up[1] - um[1]) / (2 * h)
    sxx, syy, sxy, *_ = sigma_cart(x, y, A)
    nx = np.cos(th); ny = np.sin(th)
    tx = sxx * nx + sxy * ny
    ty = sxy * nx + syy * ny
    W = W_cart(x, y, A)
    integ = (W * np.cos(th) - (tx * dux_dx + ty * duy_dx)) * R
    return np.trapezoid(integ, th)


def U_numeric(R, A, nr=7, nth=40001):
    rr = np.linspace(1e-9 * R, R, nr)
    th = np.linspace(0.0, 2 * np.pi, nth)
    RR, TT = np.meshgrid(rr, th, indexing='ij')
    x = RR * np.cos(TT); y = RR * np.sin(TT)
    se = se_cart(x, y, A)
    W = (n / (n + 1.0)) * alpha * eps0 * se ** (n + 1.0) / sig0 ** n
    # U = int_0^R int_{-pi}^{pi} W r dth dr  -> integrate over full circle
    inner = np.trapezoid(W * RR, th, axis=1)     # per r
    return np.trapezoid(inner, rr)


A = 1.0            # amplitude in metres^(2-s); ratios are A-independent
R = 1e-3
Jn = J_numeric(A, R)
Un = U_numeric(R, A)
print(f"A={A}, R={R}: J_numeric={Jn:.8e}, U_numeric={Un:.8e}")
print(f"U/(J R) = {Un/(Jn*R):.10f}")

# R-independence checks
for RR in [3e-4, 1e-3, 3e-3, 1e-2]:
    print(f"  R={RR:.0e}: J={J_numeric(A,RR):.8e}  U/(JR)={U_numeric(RR,A)/(J_numeric(A,RR)*RR):.8f}")

# =====================================================================
# H : sigma_tt(r,0) = H sig0 [J/(alpha eps0 sig0 r)]^p
#   sigma_tt(r,0) = sig0 A st~(0) r^{-p}
#   => H = A st~(0) [alpha eps0 sig0 / J]^p    (r cancels; A cancels too)
# =====================================================================
sr0, st0, srt0, se0, *_ = angular(np.array([0.0]))
st0 = st0[0]
A_ref = 1.0
J_ref = J_numeric(A_ref, 1e-3)
H = A_ref * st0 * (alpha * eps0 * sig0 / J_ref) ** p
print(f"st~(0) = {st0:.12f}")
print(f"H = A*st0*(alpha eps0 sig0/J)^p = {H:.10f}")

# check r- and A-independence of H
for AR in [0.5, 1.0, 2.0]:
    for RR in [1e-3, 1e-2]:
        Jx = J_numeric(AR, RR)
        print(f"  A={AR} R={RR}: H={AR*st0*(alpha*eps0*sig0/Jx)**p:.10f}")

# =====================================================================
# CONSISTENCY CHECK: does the dimensional field satisfy the constitutive
# law  eps_e = alpha eps0 (se/sig0)^n ?   Compute the Cartesian strain
# tensor from u_cart by finite differences, form eps_e, and compare with
# alpha eps0 (se/sig0)^n using se from the stress field.  If the field is
# self-consistent, the ratio must be 1 everywhere.
# =====================================================================
def strain_cart(x, y, A, h=1e-9):
    def u(xx, yy):
        return u_cart(xx, yy, A)
    up = u(x + h, y); um = u(x - h, y)
    vp = u(x, y + h); vm = u(x, y - h)
    uxx = (up[0] - um[0]) / (2 * h)
    uyy = (vp[1] - vm[1]) / (2 * h)
    uxy = 0.5 * ((up[1] - um[1]) / (2 * h) + (vp[0] - vm[0]) / (2 * h))
    return uxx, uyy, uxy

def eps_e_from_u(x, y, A):
    exx, eyy, exy = strain_cart(x, y, A)
    # plane strain: e_zz = 0, deviatoric strain = sqrt( (2/3) e'_ij e'_ij )
    em = (exx + eyy) / 3.0
    dxx, dyy, dxy = exx - em, eyy - em, -em
    return np.sqrt((2.0 / 3.0) * (dxx ** 2 + dyy ** 2 + em ** 2 + 2 * dxy ** 2))

def eps_e_from_stress(x, y, A):
    se = se_cart(x, y, A)
    return alpha * eps0 * (se / sig0) ** n

print("\nConstitutive consistency eps_e(u) vs alpha eps0 (se/sig0)^n:")
for (xx, yy) in [(1e-3, 1e-4), (5e-4, 3e-4), (2e-3, 1e-3), (1e-3, 5e-4)]:
    e_u = eps_e_from_u(np.array([xx]), np.array([yy]), A)[0]
    e_s = eps_e_from_stress(np.array([xx]), np.array([yy]), A)[0]
    print(f"  (x,y)=({xx:.0e},{yy:.0e}): eps_e(u)={e_u:.6e}  eps_e(stress)={e_s:.6e}  ratio={e_u/e_s:.6f}")

# =====================================================================
# FINAL: p, H, C
# =====================================================================
C = Un / (Jn * R)
# =====================================================================
# Closed-form check: eq.(24) I of the SAME angular field, and the
# resulting closed-form J = alpha eps0/sig0^n A^(n+1) I.  Compare to the
# validated numeric J.
# =====================================================================
th = np.linspace(0.0, np.pi, 400001)
sr, st, srt, se, ur, urp, ut, utp = angular(th)
integrand = ((n / (n + 1.0)) * se ** (n + 1.0) * np.cos(th)
             - (np.sin(th) * (sr * (ut - urp) - srt * (ur + utp))
                + aexp * np.cos(th) * (sr * ur + srt * ut)))
I_eq24 = 2.0 * np.trapezoid(integrand, th)
J_closed = alpha * eps0 / sig0 ** n * A ** (n + 1.0) * I_eq24
print("\n==== CLOSED-FORM vs NUMERIC J ====")
print(f"eq.(24) I (this field)   = {I_eq24:.8f}")
print(f"J_closed = a e0/s0^n A^(n+1) I = {J_closed:.6e}")
print(f"J_numeric                      = {Jn:.6e}")
print(f"ratio J_closed/J_numeric       = {J_closed/Jn:.8f}")

print("\n==== FINAL RESULTS ====")
print(f"p = {p:.12f}")
print(f"H = {H:.10f}")
print(f"C = {Un/(Jn*R):.10f}")

# =====================================================================
# Independent verification 1: the r^{-p} law and H from sigma_tt(r,0).
#   sigma_tt(r,0) should equal H sig0 [J/(alpha eps0 sig0 r)]^p for ANY r,
#   using the SAME J.  Check at three radii.
# =====================================================================
print("\n[verify] sigma_tt(r,0) vs H sig0 [J/(alpha eps0 sig0 r)]^p:")
for RR in [1e-4, 1e-3, 1e-2]:
    sth = sigma_cart(np.array([RR]), np.array([0.0]), A)
    stt_val = sth[4][0]                       # sigma_theta_theta = ST at theta=0
    pred = H * sig0 * (J_ref / (alpha * eps0 * sig0 * RR)) ** p
    print(f"  r={RR:.0e}: sigma_tt={stt_val:.6e}  H-formula={pred:.6e}  ratio={stt_val/pred:.8f}")

# Independent verification 2: check the exponent directly from two radii
r1, r2 = 1e-4, 1e-3
s1 = sigma_cart(np.array([r1]), np.array([0.0]), A)[4][0]
s2 = sigma_cart(np.array([r2]), np.array([0.0]), A)[4][0]
pexp = -np.log(s2 / s1) / np.log(r2 / r1)
print(f"\n[verify] measured exponent from sigma_tt(r,0): p = {pexp:.12f} (expect {p:.12f})")

# Independent verification 3: U(R) from the closed form
_tg = np.linspace(0, np.pi, 40001)
shape2 = 2.0 * np.trapezoid(angular(_tg)[3] ** (n + 1.0), _tg)
U_closed = (n / (n + 1.0)) * alpha * eps0 / sig0 ** n * A ** (n + 1.0) * R * shape2
print(f"\n[verify] U numeric={Un:.6e}  U closed-form={U_closed:.6e}  ratio={Un/U_closed:.8f}")