import numpy as np
import sys
sys.path.insert(0, "knowledge/checks")
from crack_tip_checks import check_field, format_report

E = 70000.0; nu = 0.30; sy = 10.0; J = 0.800; Czz0 = -0.005
A_hyd = (1 - 2 * nu) / (3 * E)
def A_dev(T): return 1.5 / T + (nu - 0.5) / E
def eff_props(T):
    G = 1.0 / (2.0 * A_dev(T)); K = 1.0 / (9.0 * A_hyd)
    Eeff = 9 * K * G / (3 * K + G); nueff = (3 * K - 2 * G) / (2 * (3 * K + G))
    return G, K, Eeff, nueff
from scipy.optimize import brentq
def Czz_of_T(T):
    G, K, Eeff, nueff = eff_props(T)
    KA = np.sqrt(J * Eeff / (2 * np.pi)); return -2 * nueff / Eeff * KA
T = brentq(lambda T: Czz_of_T(T) - Czz0, 1000, 60000, xtol=1e-13, rtol=1e-15)
G, K, Eeff, nueff = eff_props(T)
KA = np.sqrt(J * Eeff / (2 * np.pi))
print("T=%.9f Eeff=%.9f nueff=%.9f KA=%.9f" % (T, Eeff, nueff, KA))

# ---------------- check_field (linear elastic effective material) ----------------
tt = np.linspace(0, np.pi, 4001)
c2, c32 = np.cos(tt / 2), np.cos(1.5 * tt)
s2, s32 = np.sin(tt / 2), np.sin(1.5 * tt)
Srr = 1.25 * c2 - 0.25 * c32
Stt = 0.75 * c2 + 0.25 * c32
Srt = 0.25 * s2 + 0.25 * s32
rep = check_field(tt, Srr, Stt, Srt, p=0.5, n=1.0, beta=np.pi, plane="stress",
                  law="linear", nu=nueff, mode="I", leading=True, tol=1e-3, tol_eq=1e-5)
print("\n=== check_field (effective linear, nu_eff=%.6f) ===" % nueff)
print(format_report(rep))

# ---------------- direct J integral on circles ----------------
kappa = (3 - nueff) / (1 + nueff)
def u_cart(x, y):
    r = np.hypot(x, y); t = np.arctan2(y, x)
    ur = (KA / (2 * G)) * np.sqrt(r) * (kappa - np.cos(t)) * np.cos(t / 2)
    ut = -(KA / (2 * G)) * np.sqrt(r) * (kappa - np.cos(t)) * np.sin(t / 2)
    return ur * np.cos(t) - ut * np.sin(t), ur * np.sin(t) + ut * np.cos(t)
def stress_cart(x, y):
    r = np.hypot(x, y); t = np.arctan2(y, x)
    srr = KA * (1.25 * np.cos(t/2) - 0.25 * np.cos(1.5*t))
    stt = KA * (0.75 * np.cos(t/2) + 0.25 * np.cos(1.5*t))
    srt = KA * (0.25 * np.sin(t/2) + 0.25 * np.sin(1.5*t))
    c, s = np.cos(t), np.sin(t)
    sxx = srr*c**2 + stt*s**2 - 2*srt*s*c
    syy = srr*s**2 + stt*c**2 + 2*srt*s*c
    sxy = (srr-stt)*s*c + srt*(c**2-s**2)
    return sxx/np.sqrt(r), syy/np.sqrt(r), sxy/np.sqrt(r)
def strain_cart(x, y):
    r = np.hypot(x, y); t = np.arctan2(y, x)
    sxx, syy = None, None
    srr = KA*(1.25*np.cos(t/2)-0.25*np.cos(1.5*t)); stt = KA*(0.75*np.cos(t/2)+0.25*np.cos(1.5*t))
    srt = KA*(0.25*np.sin(t/2)+0.25*np.sin(1.5*t))
    exx = (srr - nueff*stt)/Eeff
    eyy = (stt - nueff*srr)/Eeff
    exy = (1+nueff)*srt/Eeff
    # rotate polar-axes strains to cartesian
    c, s = np.cos(t), np.sin(t)
    X = exx; Y = eyy; XY = exy  # these are in polar axes; convert
    exx_c = X*c**2 + Y*s**2 - 2*XY*s*c
    eyy_c = X*s**2 + Y*c**2 + 2*XY*s*c
    exy_c = (X-Y)*s*c + XY*(c**2-s**2)
    return exx_c/np.sqrt(r), eyy_c/np.sqrt(r), exy_c/np.sqrt(r)

def J_circle(R, k=200001):
    a = np.linspace(-np.pi+1e-6, np.pi-1e-6, k)
    x, y = R*np.cos(a), R*np.sin(a)
    sxx, syy, sxy = stress_cart(x, y)
    exx, eyy, exy = strain_cart(x, y)
    W = 0.5*(sxx*exx + syy*eyy + 2*sxy*exy)
    h = 1e-6*R
    u1 = np.array(u_cart(x+h, y)); u0 = np.array(u_cart(x-h, y))
    dux = (u1[0]-u0[0])/(2*h); duy = (u1[1]-u0[1])/(2*h)
    Nx, Ny = np.cos(a), np.sin(a)
    tx, ty = sxx*Nx + sxy*Ny, sxy*Nx + syy*Ny
    return np.trapezoid((W*Nx - (tx*dux + ty*duy)) * R, a)

print("\n=== direct J integral ===")
for R in (0.05, 0.1, 0.2, 1.0):
    print("  R=%.3f  J=%.10f" % (R, J_circle(R)))
print("  closed form J = K_A^2*2*pi/Eeff = %.10f" % (2*np.pi*KA**2/Eeff))

# ---------------- part 3 : C_y(90) ----------------
t90 = np.pi/2
Ur90 = (KA/(2*G))*(kappa - np.cos(t90))*np.cos(t90/2)
print("\nPart 3: C_y(90) = %.9f mm^0.5" % Ur90)

# independent: integrate du_r/dr = eps_rr along theta=90 from 0
def err_ang(t):
    srr = KA*(1.25*np.cos(t/2)-0.25*np.cos(1.5*t)); stt = KA*(0.75*np.cos(t/2)+0.25*np.cos(1.5*t))
    return (srr - nueff*stt)/Eeff
print("  via 2*err(90) = %.9f" % (2*err_ang(t90)))

# ---------------- part 4 : screening interval on theta=0 ----------------
print("\nPart 4: screens on theta=0")
sig0 = lambda r: KA/np.sqrt(r)          # sig_yy = sig_e at theta=0
ra = (KA/(10*sy))**2                     # sig_e >= 10 sy
exx0 = lambda r: (1-nueff)/Eeff * sig0(r)
ezz0 = lambda r: -2*nueff/Eeff * sig0(r)
M = max(abs((1-nueff)/Eeff), abs(2*nueff/Eeff))
rb = (KA*M/0.05)**2                      # max principal strain <= 0.05
print("  principal strains at theta=0: exx=%.6e*r^-1/2, ezz=%.6e*r^-1/2" % ((1-nueff)/Eeff*KA, -2*nueff/Eeff*KA))
print("  r_a (sig_e>=10 sy)       = %.9f mm" % ra)
print("  r_b (max|eps_princ|<=5%%) = %.9f mm" % rb)
print("  r_c (r<=0.2)             = 0.2 mm")
lo, hi = rb, min(ra, 0.2)
print("  interval = [%.9f, %.9f] mm  nonempty=%s" % (lo, hi, lo <= hi))
print("  checks: sig_e(lo)=%.6f  sig_e(hi)=%.6f  maxstrain(lo)=%.6f maxstrain(hi)=%.6f"
      % (sig0(lo), sig0(hi), M*sig0(lo), M*sig0(hi)))
