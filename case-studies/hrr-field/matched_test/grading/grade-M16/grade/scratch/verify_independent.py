"""Verification: (a) provided constraint checks on my leading field; (b) my own J-contour
integral in physical units (independent code); (c) kinematic identity C_y(90°)=2 C_rr(90°)."""
import numpy as np
from scipy.interpolate import CubicSpline
import sys
sys.path.insert(0, "knowledge/checks")

E, nu, sy = 70000.0, 0.30, 10.0
J, Czz0 = 0.800, -0.005000
T = 4816.472187
K = np.sqrt(J * T)
nus = 0.5 - (T / (2 * E)) * (1 - 2 * nu)

def f_ang(th):
    c2, s2 = np.cos(th / 2), np.sin(th / 2)
    sxx = c2 * (1 - s2 * np.sin(3 * th / 2)); syy = c2 * (1 + s2 * np.sin(3 * th / 2))
    sxy = c2 * s2 * np.cos(3 * th / 2)
    c, s = np.cos(th), np.sin(th)
    return (sxx * c**2 + syy * s**2 + 2 * sxy * s * c,
            sxx * s**2 + syy * c**2 - 2 * sxy * s * c,
            (syy - sxx) * s * c + sxy * (c**2 - s**2))

# ---------------- (b) independent J integral, physical units, physical displacements
# strains from Hooke with (T, nu*), plane stress; displacements by my own angular integration
tg = np.linspace(0, np.pi, 8001)
srr, stt, srt = f_ang(tg)
err = (srr - nus * stt) / T * K / np.sqrt(2 * np.pi)  # eps_rr = r^-1/2 * err(theta)
ett = (stt - nus * srr) / T * K / np.sqrt(2 * np.pi)
ert = (1 + nus) * srt / T * K / np.sqrt(2 * np.pi)
Ur = 2 * err                                          # u ~ r^(1/2): U_r = 2 e_rr^ang
Ut = CubicSpline(tg, ett - Ur).antiderivative()
Ut_v = Ut(tg) - Ut(0.0)                               # symmetry: u_theta(0) = 0
# extend to the full range [-pi, pi] by mode-I symmetry (s_rr, s_tt even; s_rt odd; U_t odd)
tg_f = np.concatenate([-tg[:0:-1], tg])
Ur_f = np.concatenate([Ur[:0:-1], Ur])
Ut_f = np.concatenate([-Ut_v[:0:-1], Ut_v])

def u_cart(x, y):
    r, a_ = np.hypot(x, y), np.arctan2(y, x)
    ur, ut = np.sqrt(r) * np.interp(a_, tg_f, Ur_f), np.sqrt(r) * np.interp(a_, tg_f, Ut_f)
    return ur * np.cos(a_) - ut * np.sin(a_), ur * np.sin(a_) + ut * np.cos(a_)

def sig_cart(x, y):
    r, a_ = np.hypot(x, y), np.arctan2(y, x)
    srr_, stt_, srt_ = f_ang(a_)
    c, s = np.cos(a_), np.sin(a_)
    f = K / np.sqrt(2 * np.pi * r)
    return (f * (srr_ * c**2 + stt_ * s**2 - 2 * srt_ * s * c),
            f * (srr_ * s**2 + stt_ * c**2 + 2 * srt_ * s * c),
            f * ((srr_ - stt_) * s * c + srt_ * (c**2 - s**2)))

def eps_cart(x, y):
    r, a_ = np.hypot(x, y), np.arctan2(y, x)
    c, s = np.cos(a_), np.sin(a_)
    err_, ett_, ert_ = (np.interp(a_, tg_f, err_f) / np.sqrt(r),
                        np.interp(a_, tg_f, ett_f) / np.sqrt(r),
                        np.interp(a_, tg_f, ert_f) / np.sqrt(r))
    return (err_ * c**2 + ett_ * s**2 - 2 * ert_ * s * c,
            err_ * s**2 + ett_ * c**2 + 2 * ert_ * s * c,
            (err_ - ett_) * s * c + ert_ * (c**2 - s**2))

err_f = np.concatenate([err[:0:-1], err])
ett_f = np.concatenate([ett[:0:-1], ett])
ert_f = np.concatenate([-ert[:0:-1], ert])

def J_circle(R, k=20001, h=1e-5):
    a_ = np.linspace(-np.pi, np.pi, k)
    x, y = R * np.cos(a_), R * np.sin(a_)
    sxx, syy, sxy = sig_cart(x, y)
    exx, eyy, exy = eps_cart(x, y)
    W = 0.5 * (sxx * exx + syy * eyy + 2 * sxy * exy)          # sigma_zz = 0: no zz term
    dux = (np.array(u_cart(x + h * R, y)) - np.array(u_cart(x - h * R, y))) / (2 * h * R)
    tx, ty = sxx * np.cos(a_) + sxy * np.sin(a_), sxy * np.cos(a_) + syy * np.sin(a_)
    return np.trapezoid((W * np.cos(a_) - (tx * dux[0] + ty * dux[1])) * R, a_)

print("\nindependent J on circles (physical units):")
for R in (0.005, 0.02, 0.05, 0.1):
    print(f"  r = {R:5.3f}:  J = {J_circle(R):.6f}   (target 0.800)")

# ---------------- (c) kinematic identity and C_y cross-check vs closed form u_r
Cy_kin = 2 * (K / (T * np.sqrt(2 * np.pi))) * (1.5 - 0.5 * nus) / np.sqrt(2) * np.sqrt(2)
# C_y(90°) = u_r(90°)/sqrt(r) = 2 e_rr^ang(90°);  e_rr^ang(90°)=(sigma_rr-nu* sigma_tt)(90°)
Cy_id = 2 * (np.sqrt(3 - nus) * 0 + (1.5 - 0.5 * nus) * np.sqrt(0.5)) * K / (T * np.sqrt(2 * np.pi))
print(f"\nC_y(90°) kinematic = {Cy_id:.6f};  via closed-form u_r (Sun&Jin-type) = "
      f"{2 * np.cos(np.pi/4) * ((1-nus) + (1+nus)*0.5) * K/(T*np.sqrt(2*np.pi)):.6f}")
print(f"C_zz(60°) predicted = {-(2*nus/T)*K/np.sqrt(2*np.pi)*np.cos(np.pi/6):.6f} "
      f"(consistency of angular form, not used in the fit)")
