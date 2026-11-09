import numpy as np
from scipy.optimize import brentq

E = 70000.0
nu = 0.30
sy = 10.0
J = 0.800
Czz_meas = -0.005000

def nustar(T):
    return 0.5 + (nu - 0.5) * T / E

# Near-tip effective (fully-plastic tangent) constitutive law:
#   E* = T ,  nu* = 1/2 + (nu-1/2) T/E
# Mode I plane-stress LEFM field with moduli (T, nu*):
#   K = sqrt(J*T)   (J = K^2/E', E'=T for plane stress)
#   Czz(0) = -(nu*/T)(sxx+syy)*sqrt(r) = -2 nu* sqrt(J/(2 pi T))
def f_T(T):
    return 2.0 * nustar(T) * np.sqrt(J / (2.0 * np.pi * T)) - abs(Czz_meas)

T = brentq(f_T, 1e-6, E * (1 - 1e-12), xtol=1e-12, rtol=1e-14)
nus = nustar(T)
K = np.sqrt(J * T)
print("T   = %.6f MPa" % T)
print("nu* = %.8f" % nus)
print("E*  = T = %.6f MPa" % T)
print("K   = %.8f N/mm^1.5" % K)

# ---- mode I angular functions (Cartesian), sigma = K/sqrt(2 pi r) * S_ij(theta)
def S_cart(th):
    c, s = np.cos(th), np.sin(th)
    ch, sh = np.cos(th/2), np.sin(th/2)
    s32 = np.sin(1.5*th)
    Sxx = ch * (1 - sh * s32)
    Syy = ch * (1 + sh * s32)
    Sxy = ch * sh * np.cos(1.5*th)
    return Sxx, Syy, Sxy

def S_polar(th):
    Sxx, Syy, Sxy = S_cart(th)
    c, s = np.cos(th), np.sin(th)
    Srr = Sxx*c*c + Syy*s*s + 2*Sxy*s*c
    Stt = Sxx*s*s + Syy*c*c - 2*Sxy*s*c
    Srt = (Syy - Sxx)*s*c + Sxy*(c*c - s*s)
    return Srr, Stt, Srt

pref = K / np.sqrt(2.0*np.pi)   # sqrt(r)*sigma = pref * S_ij

th = np.deg2rad(60.0)
Srr, Stt, Srt = S_polar(th)
Crr = pref / T * (Srr - nus*Stt)
Ctt = pref / T * (Stt - nus*Srr)
Cgr = pref / T * 2.0*(1.0+nus) * Srt
print("theta=60: Crr = %.8f, Ctt = %.8f, C_grthet = %.8f mm^0.5" % (Crr, Ctt, Cgr))

# ---- C_y(90): mode I plane stress displacement
# u_y = (K/2G) sqrt(r/2pi) sin(th/2)[kappa+1-2cos^2(th/2)], kappa=(3-nu)/(1+nu)
# C_y = K(3-nu*)/(2 T sqrt(pi))
Cy = K*(3.0-nus)/(2.0*T*np.sqrt(np.pi))
print("C_y(90) = %.8f mm^0.5" % Cy)

# cross-check C_y from finite-difference of the displacement field
G = T/(2*(1+nus))
kap = (3-nus)/(1+nus)
def uy(r, th):
    return K/(2*G)*np.sqrt(r/(2*np.pi))*np.sin(th/2)*(kap+1-2*np.cos(th/2)**2)
r0 = 1e-4
print("C_y FD check:", uy(r0, np.deg2rad(90))/np.sqrt(r0))

# ---- verification: numerical J integral on a circular arc, leading elastic field
def stress_polar(r, th):
    Srr, Stt, Srt = S_polar(th)
    fac = K/np.sqrt(2*np.pi*r)
    return fac*Srr, fac*Stt, fac*Srt

def eps_polar(r, th):
    Srr, Stt, Srt = stress_polar(r, th)
    err = (Srr - nus*Stt)/T
    ett = (Stt - nus*Srr)/T
    ert = (1+nus)/T*Srt
    return err, ett, ert

def disp_polar(r, th):
    # mode I cartesian displacements -> polar components
    ur_c = K/(2*G)*np.sqrt(r/(2*np.pi))*np.cos(th/2)*(kap-1+2*np.sin(th/2)**2)
    ut_c = K/(2*G)*np.sqrt(r/(2*np.pi))*np.sin(th/2)*(kap+1-2*np.cos(th/2)**2)
    # (ur_c, ut_c) already are u_x,u_y; rotate to polar basis
    c, s = np.cos(th), np.sin(th)
    ur = ur_c*c + ut_c*s
    ut = -ur_c*s + ut_c*c
    return ur, ut

def Jint(r, n=200000):
    ths = np.linspace(-np.pi, np.pi, n)
    fr = np.zeros(n)
    W = np.zeros(n)
    for i, th in enumerate(ths):
        Srr, Stt, Srt = stress_polar(r, th)
        err, ett, ert = eps_polar(r, th)
        srr, stt, srt = stress_polar(r+1e-7*r, th)  # not used
        # W = 1/2 sigma:eps (effective linear)
        W[i] = 0.5*(Srr*err + Stt*ett + 2.0*Srt*ert)
        # displacement gradients w.r.t. x
        h = 1e-6*r
        uxp = disp_polar(r, th)  # placeholder
    # compute u_{i,x} via finite differences in x
    # x = r cos th ; hold th fixed
    fr = np.zeros(n)
    for i, th in enumerate(ths):
        c, s = np.cos(th), np.sin(th)
        # du/dx at fixed y: use cartesian parametrization; easier: numerical gradient in x,y
        # point
        xp, yp = r*c, r*s
        h = 1e-7*r
        def ucart(x, y):
            rr = np.hypot(x, y); tt = np.arctan2(y, x)
            if rr == 0: rr = 1e-30
            ur, ut = disp_polar(rr, tt)
            cc, ss = np.cos(tt), np.sin(tt)
            ux = ur*cc - ut*ss
            uy = ur*ss + ut*cc
            return ux, uy
        ux1, uy1 = ucart(xp+h, yp)
        ux0, uy0 = ucart(xp-h, yp)
        ux_x = (ux1-ux0)/(2*h)
        uy_x = (uy1-uy0)/(2*h)
        ur, ut = disp_polar(r, th)
        Nx = np.cos(th)
        tx = Srr*Nx + Srt*s
        ty = Srt*Nx + Stt*s
        fr[i] = W[i]*Nx - (tx*ux_x + ty*uy_x)
    # integrate over ds = r dth
    return np.trapezoid(fr, ths)*r

for rtest in [0.01, 0.05, 0.15]:
    print("J numeric at r=%.3f mm: %.6f" % (rtest, Jint(rtest, 20001)))

print("---- consistency check: strain from stress-law vs finite diff of u ----")
def ucart(x, y):
    rr = np.hypot(x, y); tt = np.arctan2(y, x)
    ur, ut = disp_polar(rr, tt)
    cc, ss = np.cos(tt), np.sin(tt)
    return ur*cc - ut*ss, ur*ss + ut*cc

def eps_from_u(r, th):
    x, y = r*np.cos(th), r*np.sin(th)
    h = 1e-6*r
    ux_xp, uy_xp = ucart(x+h, y); ux_xm, uy_xm = ucart(x-h, y)
    ux_yp, uy_yp = ucart(x, y+h); ux_ym, uy_ym = ucart(x, y-h)
    uxx=(ux_xp-ux_xm)/(2*h); uxy=(ux_yp-ux_ym)/(2*h)
    uyx=(uy_xp-uy_xm)/(2*h); uyy=(uy_yp-uy_ym)/(2*h)
    exx=uxx; eyy=uyy; exy=0.5*(uxy+uyx)
    c,s=np.cos(th),np.sin(th)
    err=exx*c*c+eyy*s*s+2*exy*s*c
    ett=exx*s*s+eyy*c*c-2*exy*s*c
    ert=(eyy-exx)*s*c+exy*(c*c-s*s)
    return err,ett,ert

for thd in [10,60,90,150]:
    th=np.deg2rad(thd)
    err,ett,ert=eps_polar(0.1,th)
    erru,ettu,ertu=eps_from_u(0.1,th)
    print(thd, "stress-law",err,ett,ert, " FD-u",erru,ettu,ertu)
