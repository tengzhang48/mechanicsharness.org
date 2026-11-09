"""Infer T and leading-field coefficients for the slit-tip problem (TASK.md).

Units: mm, N, MPa. Plane stress (sigma_zz = 0), eps_zz free.
Total-strain law: eps = (1+nu)/E sig - nu/E tr(sig) I + (3/2)(1/T-1/E)<1-sy/se> s.
As sigma_e -> inf (r->0), <1-sy/se> -> 1, so
    eps_ij = a sig_ij + b tr(sig) d_ij,
    a = (3/2)/T + (nu-1/2)/E,  b = (1/2-nu)/E - 1/(2T)
  == Hooke plane stress with E* = T, nu* = 1/2 - (1/2-nu) T/E.
Leading tip field: Williams mode I, plane stress, moduli (T, nu*), K* = sqrt(J T)
(J = K*^2/T for plane stress; J path-independent for this Hencky material).
"""
import numpy as np
from scipy.optimize import brentq
import sys, os

E, nu, sy = 70000.0, 0.30, 10.0
J = 0.800
C_zz0 = -0.005000

def ab(T):
    a = 1.5 / T + (nu - 0.5) / E
    b = (0.5 - nu) / E - 0.5 / T
    return a, b

def nustar(T):
    return 0.5 - (0.5 - nu) * T / E

# --- stress shape functions, mode I plane stress: sig_ij = K/ sqrt(2 pi r) f_ij(theta)
def f_shapes(th):
    c, s = np.cos(th / 2), np.sin(th / 2)
    s3, c3 = np.sin(3 * th / 2), np.cos(3 * th / 2)
    fxx = c * (1 - s * s3)
    fyy = c * (1 + s * s3)
    fxy = c * s * c3
    return fxx, fyy, fxy

# --- 1) solve for T from C_zz(0): tr(sig)= 2 K*/sqrt(2 pi r) at theta=0 -> eps_zz = b*2K*/sqrt(2 pi r)
def Czz0_of_T(T):
    a, b = ab(T)
    K = np.sqrt(J * T)
    return 2.0 * b * K / np.sqrt(2 * np.pi)

T = brentq(lambda T: Czz0_of_T(T) - C_zz0, 1.0, E - 1.0, xtol=1e-10, rtol=1e-14)
a, b = ab(T)
nus = nustar(T)
K = np.sqrt(J * T)

print(f"T = {T:.6f} MPa   (nu* = {nus:.6f}, E* = T, K* = {K:.6f} MPa*mm^0.5)")
print(f"check C_zz(0): {Czz0_of_T(T):.8f}  (target {C_zz0})")
print(f"linear-law identity: a-(1+nu*)/T = {a - (1 + nus) / T:.3e},  b+nu*/T = {b + nus / T:.3e}")

# --- 2) C_rr, C_tt, C_gamma_rt at theta = 60 deg
def strain_coeffs(th, T):
    """C_ij(theta) with eps_ij = C_ij / sqrt(r); returns cart components."""
    fxx, fyy, fxy = f_shapes(th)
    K_ = np.sqrt(J * T)
    pref = K_ / np.sqrt(2 * np.pi)
    Exx = pref * (a * fxx + b * (fxx + fyy))
    Eyy = pref * (a * fyy + b * (fxx + fyy))
    Exy = pref * a * fxy          # b-term has no shear part
    return Exx, Eyy, Exy

th = np.deg2rad(60.0)
Exx, Eyy, Exy = strain_coeffs(th, T)
c, s = np.cos(th), np.sin(th)
C_rr = c**2 * Exx + s**2 * Eyy + 2 * s * c * Exy
C_tt = s**2 * Exx + c**2 * Eyy - 2 * s * c * Exy
C_grt = 2 * (s * c * (Eyy - Exx) + (c**2 - s**2) * Exy)   # engineering shear coeff
print(f"\ntheta=60deg: C_rr = {C_rr:.6f}, C_tt = {C_tt:.6f}, C_gamma_rt = {C_grt:.6f} mm^0.5")

# --- 3) C_y(90deg) from the mode I plane-stress displacement field with (T, nu*)
# u_y = K(1+nu)/E sqrt(r/2pi) sin(th/2)[kappa+1-2cos^2(th/2)], kappa=(3-nu)/(1+nu) plane stress
kappa = (3 - nus) / (1 + nus)
th9 = np.deg2rad(90.0)
Cy = (1 + nus) / T * K / np.sqrt(2 * np.pi) * np.sin(th9 / 2) * (kappa + 1 - 2 * np.cos(th9 / 2) ** 2)
print(f"\nC_y(90deg) = {Cy:.6f} mm^0.5   (kappa = {kappa:.6f})")

# cross-check by integrating strains: u_y(r,90)/sqrt(r) = u_r(90)/sqrt(r) = 2*eps_rr-coeff(90)
Exx9, Eyy9, Exy9 = strain_coeffs(th9, T)
Crr90 = Eyy9  # at theta=90: e_rr = e_yy (cart)
print(f"cross-check via u_r(90deg) = 2*C_rr(90deg): {2 * Crr90:.6f}  (should match C_y)")

# second cross-check: numerically differentiate the analytic u field -> strains vs compliance
def u_cart(x, y):
    r, th = np.hypot(x, y), np.arctan2(y, x)
    ux = (1 + nus) / T * K * np.sqrt(r / (2 * np.pi)) * np.cos(th / 2) * (kappa - 1 + 2 * np.sin(th / 2) ** 2)
    uy = (1 + nus) / T * K * np.sqrt(r / (2 * np.pi)) * np.sin(th / 2) * (kappa + 1 - 2 * np.cos(th / 2) ** 2)
    return ux, uy
h = 1e-6
r0, th0 = 1.0, np.deg2rad(60)
x0, y0 = r0 * np.cos(th0), r0 * np.sin(th0)
eps_num = np.zeros((2, 2))
eps_num[0, 0] = (u_cart(x0 + h, y0)[0] - u_cart(x0 - h, y0)[0]) / (2 * h)
eps_num[1, 1] = (u_cart(x0, y0 + h)[1] - u_cart(x0, y0 - h)[1]) / (2 * h)
eps_num[0, 1] = 0.5 * ((u_cart(x0 + h, y0)[1] - u_cart(x0 - h, y0)[1]) / (2 * h)
                       + (u_cart(x0, y0 + h)[0] - u_cart(x0, y0 - h)[0]) / (2 * h))
Ex, Ey, Exy = strain_coeffs(th0, T)
eps_asym = np.array([[Ex, Exy], [Exy, Ey]]) / np.sqrt(r0)
print("displacement-gradient check at (r=1, 60deg): max diff =",
      np.max(np.abs(eps_num - eps_asym)))

# --- 4) screening interval on theta = 0
# sigma_e = K/sqrt(2 pi r)  (uniaxial-equivalent state: sig_xx=sig_yy=p, sig_zz=0 -> sigma_e = p)
# principal strains at theta=0: eps_xx=eps_yy=a*p, eps_zz=2 b p  (p = K/sqrt(2 pi r))
amp = np.array([a, a, 2 * b])  # eps_i / p
i_max = np.argmax(np.abs(amp))
slope = amp[i_max]
print(f"\nlargest |principal strain|/p = {slope:.6e} (component {['xx','yy','zz'][i_max]})")
r_max_se = K**2 / (2 * np.pi * (10 * sy) ** 2)          # sigma_e >= 10 sy
r_min_e = (K * slope / (np.sqrt(2 * np.pi) * 0.05)) ** 2  # max|eps| <= 0.05
r_lo = max(r_min_e, 0.0)
r_hi = min(r_max_se, 0.200)
print(f"sigma_e screen:  r <= {r_max_se:.6f} mm")
print(f"strain screen:   r >= {r_min_e:.6f} mm")
print(f"r <= 0.200 mm")
if r_lo <= r_hi:
    print(f"interval = [{r_lo:.6f}, {r_hi:.6f}] mm")
else:
    print("interval EMPTY")

# --- verification A: exact law at finite r -> limit of sqrt(r) eps_zz on theta=0
def eps_zz_exact(r, T):
    p = K / np.sqrt(2 * np.pi * r)                       # sig_xx = sig_yy = p
    sig = np.diag([p, p, 0.0])
    tr = 2 * p
    s = sig - tr / 3 * np.eye(3)
    se = np.sqrt(1.5 * np.einsum('ij,ij', s, s))
    mac = max(1.0 - sy / se, 0.0)
    return (1 + nu) / E * 0.0 - nu / E * tr + 1.5 * (1 / T - 1 / E) * mac * s[2, 2]

print("\nverification A: sqrt(r)*eps_zz(r,0) with the EXACT law at finite r")
for r in (1e-2, 1e-3, 1e-4, 1e-6, 1e-8):
    print(f"  r={r:9.1e}: {np.sqrt(r) * eps_zz_exact(r, T):+.7f}")

# --- verification B: eps from exact law at theta=60deg, r small -> matches C_rr etc.
def sig_field(r, th):
    fxx, fyy, fxy = f_shapes(th)
    p = K / np.sqrt(2 * np.pi * r)
    return np.array([[fxx, fxy, 0], [fxy, fyy, 0], [0, 0, 0]]) * p

def eps_exact(r, th, T):
    sig = sig_field(r, th)
    tr = sig[0, 0] + sig[1, 1]
    s = sig - tr / 3 * np.eye(3)
    se = np.sqrt(1.5 * np.einsum('ij,ij', s, s))
    mac = max(1.0 - sy / se, 0.0)
    return (1 + nu) / E * sig - nu / E * tr * np.eye(3) + 1.5 * (1 / T - 1 / E) * mac * s

r = 1e-8
th60 = np.deg2rad(60)
eps = eps_exact(r, th60, T)
c4, s4 = np.cos(th60), np.sin(th60)
Q = np.array([[c4, s4], [-s4, c4]])
ep = Q @ eps[:2, :2] @ Q.T
val = np.sqrt(r) * np.array([ep[0, 0], ep[1, 1], 2 * ep[0, 1]])
print("\nverification B: sqrt(r)*eps_polar(60deg) exact law, r=1e-8:",
      [f"{v:+.6f}" for v in val], "vs", [f"{v:+.6f}" for v in (C_rr, C_tt, C_grt)])

# --- verification C: run the provided checks on the leading field (E*=1 units, amplitude sqrt(J))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "knowledge", "checks"))
from crack_tip_checks import check_field, format_report
N = 2001
thg = np.linspace(0, np.pi, N)
Kamp = np.sqrt(J) / np.sqrt(2 * np.pi)   # so that J(E*=1) = Kamp^2 * 2 pi / ... set by S_tt(0)=K/sqrt(2pi)
fxx, fyy, fxy = f_shapes(thg)
S_rr = Kamp * (fxx * np.cos(thg)**2 + fyy * np.sin(thg)**2 + 2 * fxy * np.sin(thg) * np.cos(thg))
S_tt = Kamp * (fxx * np.sin(thg)**2 + fyy * np.cos(thg)**2 - 2 * fxy * np.sin(thg) * np.cos(thg))
S_rt = Kamp * ((fyy - fxx) * np.sin(thg) * np.cos(thg) + fxy * (np.cos(thg)**2 - np.sin(thg)**2))
# sanity: sigma_rr(0) must equal sigma_xx(0) = K/sqrt(2 pi r) coefficient, and sigma_rt(0)=0
assert abs(S_rr[0] - Kamp) < 1e-12 and abs(S_tt[0] - Kamp) < 1e-12 and abs(S_rt[0]) < 1e-12
rep = check_field(thg, S_rr, S_tt, S_rt, p=0.5, beta=np.pi, plane="stress",
                  law="linear", nu=nus, mode="I", reported_J=J, tol=1e-3, tol_eq=1e-5)
print("\nverification C: knowledge checks on leading field (linear law, nu*, plane stress)")
print(format_report(rep))

print("\n================ SUMMARY ================")
print(f"1) T                = {T:.2f} MPa")
print(f"2) C_rr(60deg)      = {C_rr:.6f} mm^0.5")
print(f"   C_tt(60deg)      = {C_tt:.6f} mm^0.5")
print(f"   C_gamma_rt(60deg)= {C_grt:.6f} mm^0.5")
print(f"3) C_y(90deg)       = {Cy:.6f} mm^0.5")
print(f"4) screen interval  = [{r_lo:.6f}, {r_hi:.6f}] mm"
      + ("" if r_lo <= r_hi else "  (EMPTY)"))
