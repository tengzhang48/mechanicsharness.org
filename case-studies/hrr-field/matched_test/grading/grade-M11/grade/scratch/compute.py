import numpy as np
from scipy.optimize import brentq

E, nu, sy = 70000.0, 0.30, 10.0
J = 0.800          # N/mm
Czz = -0.005000    # mm^1/2

K_bulk = E / (3.0 * (1.0 - 2.0 * nu))      # 3K = E/(1-2nu)
threeK = 3.0 * K_bulk

# Effective shear compliance in fully yielded state:
#   1/G_eff = 2[(1+nu)/E + 1.5(1/T - 1/E)] = 3/T + (2nu-1)/E
# Elimination of T gives G_eff = 3K*T/(3*(3K)-T) style; but solve from ratio.
# Relation from Czz and J:   E_eff/nu_eff^2 = 2J/(pi*Czz^2)
R = 2.0 * J / (np.pi * Czz**2)

def ratio(G):
    Ee = 9.0 * K_bulk * G / (threeK + G)
    nue = (threeK - 2.0 * G) / (2.0 * (threeK + G))
    return Ee / nue**2

G = brentq(lambda g: ratio(g) - R, 1.0, 80000.0, xtol=1e-10, rtol=1e-14)
Ee = 9.0 * K_bulk * G / (threeK + G)
nue = (threeK - 2.0 * G) / (2.0 * (threeK + G))
# T from 1/G = 3/T + (2nu-1)/E
T = 3.0 / (1.0 / G - (2.0 * nu - 1.0) / E)

print("G_eff  =", G)
print("E_eff  =", Ee)
print("nu_eff =", nue)
print("T      =", T)
print("E_eff - T =", Ee - T)

# Amplitude
K = np.sqrt(J * Ee)                 # K_I = sqrt(J E_eff) plane stress
A = K / np.sqrt(2.0 * np.pi)        # sqrt(r)*sigma_xx etc prefactor on theta=0
print("K_I    =", K, " A = K/sqrt(2pi) =", A)
print("Czz predicted =", -2.0 * nue / Ee * A)

kappa = (3.0 - nue) / (1.0 + nue)   # plane stress

def stress_cart(r, th, K):
    c = np.cos(th / 2.0)
    s = np.sin(th / 2.0)
    pref = K / np.sqrt(2.0 * np.pi * r)
    sxx = pref * c * (1.0 - s * np.sin(1.5 * th))
    syy = pref * c * (1.0 + s * np.sin(1.5 * th))
    sxy = pref * s * c * np.cos(1.5 * th)
    return sxx, syy, sxy

def strain_eff(sxx, syy, sxy, Ee, nue):
    exx = (sxx - nue * syy) / Ee
    eyy = (syy - nue * sxx) / Ee
    exy = (1.0 + nue) / Ee * sxy
    ezz = -nue / Ee * (sxx + syy)
    return exx, eyy, exy, ezz

th = np.deg2rad(60.0)
r = 1e-9
sxx, syy, sxy = stress_cart(r, th, K)
exx, eyy, exy, ezz = strain_eff(sxx, syy, sxy, Ee, nue)
# rotate tensor to polar
c, s = np.cos(th), np.sin(th)
Q = np.array([[c, s], [-s, c]])
sig = np.array([[sxx, sxy], [sxy, syy]])
eps = np.array([[exx, exy], [exy, eyy]])
sigp = Q @ sig @ Q.T
epsp = Q @ eps @ Q.T
srr, srt = sigp[0, 0], sigp[0, 1]
stt = sigp[1, 1]
err, ert, ett = epsp[0, 0], epsp[0, 1], epsp[1, 1]
print("\nAt theta=60 deg")
print("stress polar srr,stt,srt =", srr, stt, srt)
print("C_rr    =", np.sqrt(r) * err)
print("C_tt    =", np.sqrt(r) * ett)
print("C_grth  =", np.sqrt(r) * 2 * ert)

# closed forms
half = th / 2.0
sh = np.sin(half)**2
Crr_f = 0.0025 / nue * np.cos(half) * ((1 + sh) - nue * (1 - sh))
Ctt_f = 0.0025 / nue * np.cos(half) * ((1 - sh) - nue * (1 + sh))
Cgt_f = 0.005 * (1 + nue) / nue * np.sin(half) * np.cos(half)**2
print("closed C_rr,C_tt,C_grth =", Crr_f, Ctt_f, Cgt_f)

# C_y(90)
Cy = 0.0025 * np.sqrt(2.0) * (3.0 - nue) / (2.0 * nue)
print("C_y(90) closed =", Cy)

# direct displacement check of C_y and strains via finite difference
def disp(r, th, K):
    c, s = np.cos(th), np.sin(th)
    fac = K / (2.0 * G) * np.sqrt(r / (2.0 * np.pi))
    ux = fac * np.cos(th / 2.0) * (kappa - np.cos(th))
    uy = fac * np.sin(th / 2.0) * (kappa - np.cos(th))
    return ux, uy

# C_y(90): u_y(r,90)/sqrt(r)
r0 = 1e-9
ux, uy = disp(r0, np.pi / 2.0, K)
print("C_y(90) direct =", uy / np.sqrt(r0))

# finite-difference kinematic strains at theta=60
def u_polar(r, th):
    ux, uy = disp(r, th, K)
    ur = ux * np.cos(th) + uy * np.sin(th)
    ut = -ux * np.sin(th) + uy * np.cos(th)
    return ur, ut

r0 = 1e-6
dr = 1e-9
dth = 1e-9
ur0, ut0 = u_polar(r0, th)
urp, _ = u_polar(r0 + dr, th)
urm, _ = u_polar(r0 - dr, th)
ur_p, ut_p = u_polar(r0, th + dth)
ur_m, ut_m = u_polar(r0, th - dth)
_, utp = u_polar(r0 + dr, th)
_, utm = u_polar(r0 - dr, th)
err_fd = (urp - urm) / (2 * dr)
ett_fd = ur0 / r0 + (ut_p - ut_m) / (2 * dth) / r0
ert_fd = 0.5 * ((ur_p - ur_m) / (2 * dth) / r0 + (utp - utm) / (2 * dr) - ut0 / r0)
print("FD  C_rr,C_tt,C_grth =", np.sqrt(r0) * err_fd, np.sqrt(r0) * ett_fd, np.sqrt(r0) * 2 * ert_fd)

# ---- Part 4: screening interval on theta=0 ----
# sigma_e = sigma_xx = A/sqrt(r); principal strains: (1-nue)sig/Ee twice, -2nue sig/Ee
r1 = (A / (10.0 * sy))**2                 # sigma_e = 10 sy
M = max(1.0 - nue, 2.0 * nue)
r2 = (M * A / (0.0500 * Ee))**2           # max |principal strain| = 0.05
print("\nA =", A, " M =", M)
print("r1 (sigma_e=10sy) =", r1)
print("r2 (strain=0.05)  =", r2)
rhi = min(r1, 0.200)
print("interval =", r2, "to", rhi, " empty?", r2 > rhi)

# verification: J integral on circles for effective-linear material
def Jnum(rad, K, Ee, nue):
    n = 8000
    ths = np.linspace(-np.pi + 1e-6, np.pi - 1e-6, n)
    sxx, syy, sxy = stress_cart(rad, ths, K)
    exx, eyy, exy, ezz = strain_eff(sxx, syy, sxy, Ee, nue)
    W = 0.5 * (sxx * exx + syy * eyy + 2.0 * sxy * exy)
    c, s = np.cos(ths), np.sin(ths)
    tx = sxx * c + sxy * s
    ty = sxy * c + syy * s
    X = rad * c; Y = rad * s
    h = 1e-7 * rad
    rp = np.sqrt((X + h)**2 + Y**2); tp = np.arctan2(Y, X + h)
    rm = np.sqrt((X - h)**2 + Y**2); tm = np.arctan2(Y, X - h)
    uxp = (disp(rp, tp, K)[0] - disp(rm, tm, K)[0]) / (2 * h)
    uyp = (disp(rp, tp, K)[1] - disp(rm, tm, K)[1]) / (2 * h)
    integrand = W * c - (tx * uxp + ty * uyp)
    return np.trapezoid(integrand, ths * rad)

for rad in [1e-3, 1e-2, 0.1]:
    print("J on circle r=%.3g : %.6f" % (rad, Jnum(rad, K, Ee, nue)))
print("K^2/E_eff =", K**2 / Ee)
