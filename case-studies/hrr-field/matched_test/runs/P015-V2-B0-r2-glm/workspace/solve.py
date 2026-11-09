"""Infer tangent modulus T from C_zz(0), then compute requested coefficients.

Model: stationary mode-I crack tip, total-strain (deformation) law, plane stress
(sigma_zz = sigma_xz = sigma_yz = 0), proportional monotonic loading.
As r->0, sigma_e ~ r^-1/2 -> infinity, so <1 - sigma_y/sigma_e> -> 1 and the
leading field obeys the linear law
    eps_ij = (1+nu)/E sig_ij - nu/E sig_kk del_ij + (3/2)(1/T - 1/E) s_ij,
which is exactly Hooke's law of an isotropic material with
    E_eq = T,   nu_eq = 1/2 - (1-2 nu) T/(2E)
(bulk modulus unchanged: K = E/(3(1-2nu)); 1/E_eq = 1/(9K) + (2/3)A, A = 3/(2T)+(nu-1/2)/E).
J is path independent and, evaluated on a vanishing circle where the leading field
dominates, equals the linear-elastic value K_I^2/E_eq = K_I^2/T (plane stress).
The stress angular functions are the standard Williams mode-I functions (nu-independent).
"""
import numpy as np
from scipy.optimize import brentq
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'knowledge', 'checks'))
from crack_tip_checks import check_field, format_report

E, nu, sy = 70000.0, 0.30, 10.0
Jg, Czz_meas = 0.800, -0.005
s2p = np.sqrt(2.0 * np.pi)

# ---------------- 1. solve for T from C_zz ------------------------------------
nu_eq = lambda T: 0.5 - (1.0 - 2.0 * nu) * T / (2.0 * E)
def Czz_of_T(T):
    # eps_zz = -(nu_eq/T)(sig_xx+sig_yy); on theta=0: sig_xx+sig_yy = 2K/sqrt(2 pi r)
    K = np.sqrt(Jg * T)                      # K_I = sqrt(J * E_eq), E_eq = T
    return -2.0 * nu_eq(T) * K / (s2p * T)

fT = lambda T: Czz_of_T(T) - Czz_meas
T = brentq(fT, 10.0, 0.999 * E, xtol=1e-13, rtol=8.9e-16, maxiter=200)
nue = nu_eq(T)
K = np.sqrt(Jg * T)
fac = K / s2p / T                            # multiply angular stress functions -> strain coeff
print(f"T               = {T:.6f} MPa")
print(f"nu_eq           = {nue:.8f}")
print(f"K_I = sqrt(J T) = {K:.6f} MPa sqrt(mm)")
print(f"K/sqrt(2pi)     = {K/s2p:.6f} MPa")
print(f"check C_zz      = {Czz_of_T(T):.8f} mm^(1/2) (measured {Czz_meas})")

# ---------------- 2. Williams mode-I angular functions ------------------------
def f_cart(th):
    c2, s2 = np.cos(th / 2), np.sin(th / 2)
    s3, c3 = np.sin(3 * th / 2), np.cos(3 * th / 2)
    return c2 * (1 - s2 * s3), c2 * (1 + s2 * s3), c2 * s2 * c3
def f_polar(th):
    a, b, c = f_cart(th)
    co, si = np.cos(th), np.sin(th)
    return (a * co**2 + b * si**2 + 2 * c * si * co,
            a * si**2 + b * co**2 - 2 * c * si * co,
            (b - a) * si * co + c * (co**2 - si**2))

# closed-form polar stresses at 60 deg (literature values 5sqrt3/8, 3sqrt3/8, 3/8)
frr, ftt, frt = f_polar(np.pi / 3)
print(f"\ntheta=60deg: f_rr={frr:.8f} (5sqrt3/8={5*np.sqrt(3)/8:.8f}), "
      f"f_tt={ftt:.8f} (3sqrt3/8={3*np.sqrt(3)/8:.8f}), f_rt={frt:.8f} (3/8={3/8:.8f})")

def C_polar(th):
    frr, ftt, frt = f_polar(th)
    return ((frr - nue * ftt) * fac,
            (ftt - nue * frr) * fac,
            2.0 * (1 + nue) * frt * fac)      # engineering shear coefficient

th = np.pi / 3
crr, ctt, crt = C_polar(th)
print(f"C_rr(60deg)  = {crr:.6e} mm^(1/2)")
print(f"C_tt(60deg)  = {ctt:.6e} mm^(1/2)")
print(f"C_grt(60deg) = {crt:.6e} mm^(1/2)")

# ---------------- 3. displacement and C_y(90deg) ------------------------------
# strain field (polar): eps_rr = r^-1/2 * (f_rr - nue f_tt) * fac etc.
# integrate: u_r = r^1/2 U_r, U_r = 2*eps_rr_coeff; U_t' = eps_tt_coeff - U_r;
# 2 eps_rt coeff = U_r' - U_t/2 ; U_t(0)=0
N = 200001
tgrid = np.linspace(0.0, np.pi, N)
Frr, Ftt, Frt = f_polar(tgrid)
err = (Frr - nue * Ftt) * fac
ett = (Ftt - nue * Frr) * fac
ert = (1 + nue) * Frt * fac                    # tensor shear
Ur = 2.0 * err
dUt = ett - Ur
Ut = np.concatenate([[0.0], np.cumsum(0.5 * (dUt[1:] + dUt[:-1]) * np.diff(tgrid))])
# shear-relation residual (compatibility check of the integration)
dUr = np.gradient(Ur, tgrid)
res_shear = np.max(np.abs(dUr - 0.5 * Ut - 2.0 * ert)[N // 20:-N // 20])
print(f"\ncompatibility of integrated displacement (shear relation residual): {res_shear:.2e}")
i90 = np.searchsorted(tgrid, np.pi / 2)
Urb, Utb = np.interp(np.pi / 2, tgrid, Ur), np.interp(np.pi / 2, tgrid, Ut)
# u_y = u_r sin + u_t cos at theta=90deg
Cy_num = Urb * 1.0 + Utb * 0.0
print(f"C_y(90deg) by integrating strains from the stress field      = {Cy_num:.6e}")

# closed form (plane stress, equivalent material): u = K(1+nue)/T sqrt(r/2pi) ...
kappa = (3.0 - nue) / (1.0 + nue)
Cy_cf = K * (1 + nue) * kappa / (2.0 * np.sqrt(np.pi) * T)
print(f"C_y(90deg) closed form K(1+nu_eq)kappa/(2 sqrt(pi) E_eq)     = {Cy_cf:.6e}")

# line integral on x=0: u_y(0,y) = int eps_yy dy = 2 C_eyy(90) sqrt(y)
fxx90, fyy90, _ = f_cart(np.pi / 2)
Ceyy90 = (fyy90 - nue * fxx90) * fac
print(f"C_y(90deg) line integral 2*C_eyy(90deg)                      = {2*Ceyy90:.6e}")

# finite-difference the closed-form displacement field -> strains; compare with law
def u_closed(x, y):
    r, th = np.hypot(x, y), np.arctan2(y, x)
    c2, s2 = np.cos(th / 2), np.sin(th / 2)
    pref = K * (1 + nue) / T * np.sqrt(r) / s2p   # = K(1+nu_eq)/E_eq * sqrt(r/(2 pi))
    ux = pref * c2 * (kappa - 1 + 2 * s2**2)
    uy = pref * s2 * (kappa + 1 - 2 * c2**2)
    return ux, uy
h = 1e-6
for (rr, aa) in [(0.05, 0.6), (0.137, 2.1), (0.02, np.pi / 2)]:
    x, y = rr * np.cos(aa), rr * np.sin(aa)
    duxdx = (u_closed(x + h, y)[0] - u_closed(x - h, y)[0]) / (2 * h)
    duydy = (u_closed(x, y + h)[1] - u_closed(x, y - h)[1]) / (2 * h)
    duxdy = (u_closed(x, y + h)[0] - u_closed(x, y - h)[0]) / (2 * h)
    p, q, s = f_polar(aa)
    exx = ((p + q) - nue * (p + q)) * 0  # rebuild from cartesian stress
    a, b, c = f_cart(aa)
    sc = K / np.sqrt(2 * np.pi * rr)
    sxx, syy, sxy = a * sc, b * sc, c * sc
    exx_l = (sxx - nue * syy) / T
    eyy_l = (syy - nue * sxx) / T
    exy_l = (1 + nue) * sxy / T
    print(f"FD strains vs law at (r={rr}, th={aa:.3f}): "
          f"dUx/dx-exx={duxdx-exx_l:.2e}, dUy/dy-eyy={duydy-eyy_l:.2e}, "
          f"1/2(dUx/dy+dUy/dx)-exy={0.5*(duxdy+ (u_closed(x,y+h)[1]-u_closed(x,y-h)[1])/(2*h))-exy_l:.2e}")

# ---------------- 4. screening interval on theta=0 ----------------------------
C_in = (1.0 - nue) * fac                      # coefficient of eps_xx=eps_yy on theta=0
Czz_abs = abs(Czz_meas)
r_strain = (max(Czz_abs, C_in) / 0.05) ** 2   # max |principal strain| <= 0.05
r_sigma = K**2 / (2 * np.pi * (10 * sy)**2)   # sigma_e = K/sqrt(2 pi r) >= 10 sigma_y
r_out = 0.200
rlo, rhi = r_strain, min(r_sigma, r_out)
print(f"\nscreening: C_in={C_in:.6e}, |C_zz|={Czz_abs:.6e} -> binding {max(Czz_abs, C_in):.6e}")
print(f"r from strain screen  >= {r_strain:.6f} mm")
print(f"r from sigma_e screen <= {r_sigma:.6f} mm (K^2/(2 pi (10 sig_y)^2))")
print(f"interval on theta=0: [{rlo:.6f}, {rhi:.6f}] mm, nonempty={rlo <= rhi}")
for r in [rlo, rhi]:
    se = K / np.sqrt(2 * np.pi * r)
    ez = Czz_meas / np.sqrt(r)
    print(f"  r={r:.6f}: sigma_e={se:.3f} MPa ({se/sy:.1f} sig_y), "
          f"eps_max={max(abs(ez), C_in/np.sqrt(r)):.5f}, bracket <1-sig_y/sig_e>={1-sy/se:.4f}")

# ---------------- 5. J with the true nonlinear W on finite circles ------------
def sig_field(x, y):
    r, th = np.hypot(x, y), np.arctan2(y, x)
    a, b, c = f_cart(th)
    sc = K / np.sqrt(2 * np.pi * r)
    return a * sc, b * sc, c * sc
def W_true(sxx, syy, sxy):
    tr = sxx + syy
    We = 0.5 * ((1 + nu) / E * (sxx**2 + syy**2 + 2 * sxy**2) - nu / E * tr**2)
    sd = np.array([sxx - tr / 3, syy - tr / 3, -tr / 3, sxy])
    Se = np.sqrt(1.5 * (sd[0]**2 + sd[1]**2 + sd[2]**2 + 2 * sd[3]**2))
    Wp = np.where(Se > sy, 0.5 * (1 / T - 1 / E) * (Se**2 - sy**2), 0.0)
    return We + Wp
def J_circle(R, n=400001):
    a = np.linspace(-np.pi, np.pi, n)
    x, y = R * np.cos(a), R * np.sin(a)
    sxx, syy, sxy = sig_field(x, y)
    W = W_true(sxx, syy, sxy)
    hR = 1e-6 * R
    dux = (u_closed(x + hR, y)[0] - u_closed(x - hR, y)[0]) / (2 * hR)
    duy = (u_closed(x + hR, y)[1] - u_closed(x - hR, y)[1]) / (2 * hR)
    tx, ty = sxx * np.cos(a) + sxy * np.sin(a), sxy * np.cos(a) + syy * np.sin(a)
    return np.trapezoid((W * np.cos(a) - (tx * dux + ty * duy)) * R, a)
print("\nJ on circles (true deformation-theory W, leading-field u):")
for R in [1e-5, 1e-4, 1e-3, 6.1e-3, 0.0613, 0.2]:
    print(f"  R={R:8.5f} mm:  J = {J_circle(R):.6f} N/mm   (target 0.800)")

# ---------------- 6. full-law eps_zz check -------------------------------------
print("\nsqrt(r)*eps_zz(r,0) from the FULL law (should -> -0.005):")
for r in [1e-8, 1e-6, 1e-4, 1e-2]:
    p = K / np.sqrt(2 * np.pi * r)
    se = p
    ez = -nu / E * 2 * p + 1.5 * (1 / T - 1 / E) * max(0.0, 1 - sy / se) * (-2 * p / 3)
    print(f"  r={r:.0e}: sqrt(r)*eps_zz = {np.sqrt(r) * ez:.8f}")

# ---------------- 7. independent checks via knowledge/checks ------------------
# feed angular functions S = (K/sqrt(2pi)) f(theta) so that sigma = r^-1/2 S.
# law="linear", E=1, nu=nu_eq: strains scale by T relative to physical, so the
# field's own J equals K^2 = J*T in these units; pass reported_J = Jg*T.
tt = np.linspace(0, np.pi, 4001)
a, b, c = f_cart(tt)
co, si = np.cos(tt), np.sin(tt)
S_rr = (a * co**2 + b * si**2 + 2 * c * si * co) * K / s2p
S_tt = (a * si**2 + b * co**2 - 2 * c * si * co) * K / s2p
S_rt = ((b - a) * si * co + c * (co**2 - si**2)) * K / s2p
rep = check_field(tt, S_rr, S_tt, S_rt, p=0.5, n=1.0, beta=np.pi, plane="stress",
                  law="linear", nu=nue, mode="I", leading=True,
                  reported_J=Jg * T)
print("\n=== knowledge/checks/crack_tip_checks.py on my leading field ===")
print(format_report(rep))
print(f"\n(expected field J in E=1 units = J*T = {Jg*T:.4f})")

np.save('results.npy', np.array([T, nue, K, crr, ctt, crt, Cy_num, Cy_cf, rlo, rhi],
                                dtype=float), allow_pickle=False)
