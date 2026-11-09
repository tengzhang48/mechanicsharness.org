"""Check the angular displacement functions u~ (from integral_I) against the
plane-strain fully-plastic constitutive law  eps~e = alpha se~^n.

In the normalized field:
  sigma~ = r^(s-2) s_ang(th),  u~ = A_n r^(n(s-2)+1) u_ang(th)   (dimensionless)
  strain~  from u~ via strain-displacement (polar, incl. angular d/dth)
  eps~e = sqrt(2/3) * e'_ij e'_ij     (plane strain, e_zz=0)
Compare to alpha * se~^n.
"""
import numpy as np
from scipy.integrate import solve_ivp

n = 6.5
alpha = 0.02
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
p_shoot = -0.689256761931
sol = solve_ivp(ff, [0.0, np.pi], [1.0, 0.0, p_shoot, 0.0],
                rtol=1e-12, atol=1e-14, dense_output=True, max_step=0.002)
M = Mcoeff(s)


def ang(th):
    th = np.asarray(th, float)
    Y = sol.sol(np.atleast_1d(th.ravel()))
    se = np.sqrt(np.einsum('ij,jk,ik->i', Y[:3].T, M, Y[:3].T))
    c = 1.0 / se.max()
    Y = Y * c
    f, fp, fpp, fppp = Y
    sr = s * f + fpp; st = s * (s - 1.0) * f; srt = (1.0 - s) * fp
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
    out = [x.reshape(th.shape) for x in (sr, st, srt, se, ur, urp, ut, utp)]
    return out


# dimensionless full field, r = 1 for clean strains:
#   eps_rr = u_r,r ; eps_tt = (u_r + u_t,th)/r ; 2 eps_rt = u_r,th/r + u_t,r - u_t/r
def strain_at(th, r=1.0):
    sr, st, srt, se, ur, urp, ut, utp = ang(th)
    aexp = n * (s - 2.0) + 1.0
    ur_f = r ** aexp * ur
    ut_f = r ** aexp * ut
    err = aexp * r ** (aexp - 1.0) * ur
    ett = (ur_f + r ** aexp * utp) / r
    ert = 0.5 * (r ** aexp * urp / r + aexp * r ** (aexp - 1.0) * ut - ut_f / r)
    return err, ett, ert, se


def epse(err, ett, ert):
    em = (err + ett) / 3.0
    de = 0.0  # e_zz = 0
    return np.sqrt((2.0 / 3.0) * ((err - em) ** 2 + (ett - em) ** 2 + em ** 2
                                  + 2 * ert ** 2))


th = np.array([0.005, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
err, ett, ert, se = strain_at(th)
ee = epse(err, ett, ert)
print("theta, eps_e(u), alpha*se^n, ratio:")
for t, a, b in zip(th, ee, alpha * se ** n):
    print(f"  {t:5.3f}  {a:.6e}  {b:.6e}  {a/b:.6f}")

# =====================================================================
# Dimensionless J integral of the normalized field, using the SAME
# integrand as solve.py's eq.(24), and compare:
#   (A) the paper's eq.(24) I  (from solve.py)  ->  4.818379
#   (B) J computed with the consistent field  (u = r^aexp u_ang, so that
#       eps_e(u) = se^n exactly; the alpha sits in W only if the field is
#       alpha-normalized)
# The physical constitutive law has eps_e = alpha eps0 (se/sig0)^n, so with
# sigma = sig0 A r^(s-2) s_ang,  eps = alpha eps0 A^n r^(n(s-2)) e_ang,
# u = alpha eps0 A^n r^aexp u_ang.  Both W and t.du/dx then carry ONE alpha,
# and J = alpha*(eps0/sig0^n)*Kdim^(n+1)*Jhat with Jhat the coefficient below.
# =====================================================================
def Jhat(nth=400001):
    th = np.linspace(-np.pi, np.pi, nth)
    # build field with A=1, sigma0=1, eps0=1, alpha=1 folded into W only via
    # the ALPHA=1 dimensionless constitutive (so eps_e = se^n).
    r = 1.0
    sr, st, srt, se, ur, urp, ut, utp = ang(th)
    aexp = n * (s - 2.0) + 1.0
    # displacement components in polar, dimensionless (u = r^aexp u_ang)
    Ur = ur; Ut = ut
    # W = n/(n+1) sigma_ij eps_ij with eps_e = se^n (alpha=1)
    W = (n / (n + 1.0)) * se ** (n + 1.0)
    # traction work t_i du_i/dx on the circle of radius r=1.
    # In polar coords the derivative wrt x (=d/dx, fixed Cartesian) is
    # d/dx = cos th d/dr - (sin th /r) d/dth, and the polar basis rotates.
    # u_x = Ur cos th - Ut sin th
    ux = Ur * np.cos(th) - Ut * np.sin(th)
    uy = Ur * np.sin(th) + Ut * np.cos(th)
    # radial derivative: u = r^aexp u_ang, so du/dr = aexp * u/r
    dux_dr = aexp * ux / r
    # angular derivative
    dux_dth = (urp * np.cos(th) - ur * np.sin(th)) - (utp * np.sin(th) + ut * np.cos(th))
    duy_dth = (urp * np.sin(th) + ur * np.cos(th)) + (utp * np.cos(th) - ut * np.sin(th))
    dux_dx = np.cos(th) * dux_dr - np.sin(th) / r * dux_dth
    duy_dx = np.cos(th) * (aexp * uy / r) - np.sin(th) / r * duy_dth
    nx = np.cos(th); ny = np.sin(th)
    tx = sr * nx + srt * ny
    ty = srt * nx + st * ny
    integ = (W * np.cos(th) - (tx * dux_dx + ty * duy_dx)) * r
    return np.trapezoid(integ, th)


Jh = Jhat()
print(f"\nJhat (dimensionless J, alpha=1 field) = {Jh:.8f}")
print(f"eq.(24) I (solve.py)                  = 4.81837907")
print(f"ratio I/Jhat                          = {4.81837907/Jh:.8f}")