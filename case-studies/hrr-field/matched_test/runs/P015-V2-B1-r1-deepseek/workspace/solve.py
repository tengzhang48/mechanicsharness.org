import numpy as np
from scipy.optimize import brentq

E = 70000.0
nu = 0.30
sy = 10.0
J = 0.800
Czz = -0.005

def alpha1(T):
    return (1 + nu) / E + 1.5 * (1.0 / T - 1.0 / E)

def alpha2(T):
    return nu / E + 0.5 * (1.0 / T - 1.0 / E)

def K_of_T(T):
    return np.sqrt(J * T)

def Czz_of_T(T):
    return -2.0 * alpha2(T) * K_of_T(T) / np.sqrt(2 * np.pi)

f = lambda T: Czz_of_T(T) - Czz
T = brentq(f, 1.0, E - 1e-9, xtol=1e-10, rtol=1e-14)
K = K_of_T(T)
a1, a2 = alpha1(T), alpha2(T)

print("T   =", T)
print("K   =", K)
print("a1  =", a1, "a2 =", a2)
print("2a2 =", 2 * a2)
print("Czz(T) =", Czz_of_T(T))

# effective 2D plane-stress material
nu2 = a2 * T
mu2 = T / (2 * (1 + nu2))
kappa2 = (3 - nu2) / (1 + nu2)
print("nu2 =", nu2, "mu2 =", mu2, "kappa2 =", kappa2)
print("check E2 = 1/(a1-a2) =", 1.0 / (a1 - a2), " T =", T, " nu2 alt =", a2 / (a1 - a2))

# --- angular stress functions: sigma = (K/sqrt(2 pi r)) f(theta)
def f_cart(th):
    c, s = np.cos(th / 2), np.sin(th / 2)
    c3, s3 = np.cos(3 * th / 2), np.sin(3 * th / 2)
    fxx = c * (1 - s * s3)
    fyy = c * (1 + s * s3)
    fxy = s * c * c3
    return fxx, fyy, fxy

def f_polar(th):
    fxx, fyy, fxy = f_cart(th)
    c, s = np.cos(th), np.sin(th)
    frr = fxx * c * c + fyy * s * s + 2 * fxy * s * c
    ftt = fxx * s * s + fyy * c * c - 2 * fxy * s * c
    frt = -(fxx - fyy) * s * c + fxy * (c * c - s * s)
    return frr, ftt, frt

th = np.radians(60.0)
frr, ftt, frt = f_polar(th)
Shat = frr + ftt
pref = K / np.sqrt(2 * np.pi)
Crr = pref * (a1 * frr - a2 * Shat)
Ctt = pref * (a1 * ftt - a2 * Shat)
Cgr = pref * (2 * a1 * frt)
print("theta=60: frr,ftt,frt =", frr, ftt, frt)
print("C_rr =", Crr)
print("C_tt =", Ctt)
print("C_gr =", Cgr)

# --- C_y(90) from the standard plane-stress mode I displacement field (E2=T, nu2)
Cy = (K / (2 * mu2)) * (1 / np.sqrt(2 * np.pi)) * (np.sqrt(2) / 2) * kappa2
print("C_y(90) =", Cy)

# independent route: u_y = u_r at theta=90, u_r = 2*C_rr*r^.5
th90 = np.radians(90.0)
frr90, ftt90, frt90 = f_polar(th90)
Crr90 = pref * (a1 * frr90 - a2 * (frr90 + ftt90))
print("C_rr(90) =", Crr90, " 2*C_rr(90) =", 2 * Crr90)

# verify via direct displacement from integration of strains at theta=90 manually:
# u_r coefficient from e_rr/(1+q) with q=-1/2 => Ur = 2*e_rr*... Actually u_r = r^.5 * Ur, err = du_r/dr = .5 r^-.5 Ur
# so C_rr relates. u_theta: u_theta = r^.5 Ut, e_tt = (Ur + Ut')/... Use check below.

# --- screens on theta = 0
s0 = K / np.sqrt(2 * np.pi)  # sigma = s0 / sqrt(r)
# sigma_e = s0/sqrt(r) must be >= 10*sy
r_ey = (s0 / (10 * sy)) ** 2
# principal strains at theta=0: eps_xx=eps_yy=(a1-2a2)*sigma, eps_zz=-2a2*sigma
e_in = (a1 - 2 * a2) * s0
e_zz = -2 * a2 * s0
print("eps_in coeff", e_in, "eps_zz coeff", e_zz)
emax_coeff = max(abs(e_in), abs(e_zz))
r_strain = (emax_coeff / 0.05) ** 2
r_cap = 0.200
print("r from sigma_e>=10sy :", r_ey)
print("r from strain<=0.05  :", r_strain)
lo = r_strain
hi = min(r_ey, r_cap)
print("interval = [", lo, ",", hi, "]  empty?", hi < lo)

# --- verification: J integral on a circle using the effective plane-stress material (E2,nu2)
def u_components(r, thg):
    # standard mode I plane stress displacement
    q = np.sqrt(r)
    c, s = np.cos(thg / 2), np.sin(thg / 2)
    fac = K / (2 * mu2) / np.sqrt(2 * np.pi)
    ux = fac * q * c * (kappa2 - 1 + 2 * s * s)
    uy = fac * q * s * (kappa2 + 1 - 2 * c * c)
    return ux, uy

def stress_cart(r, thg):
    fxx, fyy, fxy = f_cart(thg)
    p = K / np.sqrt(2 * np.pi * r)
    return p * fxx, p * fyy, p * fxy

def strain_cart(r, thg):
    frr, ftt, frt = f_polar(thg)
    Shat = frr + ftt
    pref = K / np.sqrt(2 * np.pi * r)
    err = pref * (a1 * frr - a2 * Shat)
    ett = pref * (a1 * ftt - a2 * Shat)
    ert = pref * (a1 * frt)
    c, s = np.cos(thg), np.sin(thg)
    exx = err * c * c + ett * s * s - 2 * ert * s * c
    eyy = err * s * s + ett * c * c + 2 * ert * s * c
    exy = (err - ett) * s * c + ert * (c * c - s * s)
    return exx, eyy, exy

def u_cart_xy(x, y):
    r = np.hypot(x, y)
    thg = np.arctan2(y, x)
    ux, uy = u_components(r, thg)
    return ux, uy

def J_circle(R, npts=200001):
    a = np.linspace(-np.pi, np.pi, npts)
    x, y = R * np.cos(a), R * np.sin(a)
    sxx, syy, sxy = stress_cart(R, a)
    exx, eyy, exy = strain_cart(R, a)
    W = 0.5 * (sxx * exx + syy * eyy + 2 * sxy * exy)
    h = 1e-6
    uxp, uyp = u_cart_xy(x + h, y)
    uxm, uym = u_cart_xy(x - h, y)
    dux = (uxp - uxm) / (2 * h)
    duy = (uyp - uym) / (2 * h)
    Nx, Ny = np.cos(a), np.sin(a)
    tx, ty = sxx * Nx + sxy * Ny, sxy * Nx + syy * Ny
    integrand = (W * Nx - (tx * dux + ty * duy)) * R
    return np.trapezoid(integrand, a)

for R in (0.1, 0.5, 1.0, 2.0):
    print("J_contour R=", R, "=", J_circle(R), " K^2/T =", K * K / T)
