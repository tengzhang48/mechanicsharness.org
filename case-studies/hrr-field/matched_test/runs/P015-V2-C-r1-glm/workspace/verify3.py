import numpy as np
from scipy.optimize import brentq

trapz = getattr(np, "trapezoid", None) or np.trapz

E, nu, sy, J = 70000.0, 0.30, 10.0, 0.800
Czz_meas = -0.005
Tsol = brentq(lambda T: -(0.5 - (0.5 - nu) * T / E) / T * 2.0 * np.sqrt(J * T) / np.sqrt(2 * np.pi) - Czz_meas,
              100.0, 69999.0, xtol=1e-13, rtol=1e-15)
Et, nt, K = Tsol, 0.5 - (0.5 - nu) * Tsol / E, np.sqrt(J * Tsol)
phi0 = K / (Et * np.sqrt(2 * np.pi))
print(f"T={Tsol:.6f} MPa, E_t={Et:.4f}, nu_t={nt:.8f}, K_t={K:.6f}, phi0={phi0:.8e}")

# ---------------- leading field, closed forms (tangent moduli) ----------------
# stress coefficients f_ij(theta) multiply K/sqrt(2 pi r)
def ang(th):
    th = np.asarray(th, float)
    c, s = np.cos(th / 2), np.sin(th / 2)
    s3, c3 = np.sin(1.5 * th), np.cos(1.5 * th)
    fxx = c * (1 - s * s3); fyy = c * (1 + s * s3); fxy = c * s * c3
    srr, stt, srt = c * (1 + s * s), c * (1 - s * s), s * c * c
    grr = srr - nt * stt          # eps_rr coefficient / phi0
    gtt = stt - nt * srr          # eps_theta-theta coefficient / phi0
    grt = (1 + nt) * srt          # eps_rtheta (tensor) coefficient / phi0
    gxx, gyy, gxy = fxx - nt * fyy, fyy - nt * fxx, (1 + nt) * fxy
    return dict(fxx=fxx, fyy=fyy, fxy=fxy, srr=srr, stt=stt, srt=srt,
                grr=grr, gtt=gtt, grt=grt, gxx=gxx, gyy=gyy, gxy=gxy)

# displacement: u = phi0 sqrt(r) [ V_r e_r + V_th e_th ], from
#   eps_rr = dU_r/dr        => U_r = 2 phi0 g_rr
#   eps_rth = (1/2)(U_r'/r ... ) => U_th = 4 phi0 (g_rr' - g_rth)
def Vfun(th):
    th = np.asarray(th, float)
    c, s = np.cos(th / 2), np.sin(th / 2)
    grr = c * ((1 - nt) + (1 + nt) * s * s)
    dgrr = -(s / 2) * ((1 - nt) + (1 + nt) * s * s) + (1 + nt) * s * c * c
    grt = (1 + nt) * s * c * c
    return 2 * grr, 4 * (dgrr - grt), grr, dgrr, grt

def Fcart(th):
    # F = R(theta) V(theta): u = phi0 sqrt(r) F(theta) in fixed x,y basis
    ct, st = np.cos(th), np.sin(th)
    Vr, Vt, _, _, _ = Vfun(th)
    return Vr * ct - Vt * st, Vr * st + Vt * ct

def Fcart_d(th, h=1e-6):
    Fp = (Fcart(th + h)[0] - Fcart(th - h)[0], Fcart(th + h)[1] - Fcart(th - h)[1])
    Fm = (Fcart(th - h)[0], Fcart(th - h)[1])
    Pp = (Fcart(th + h)[0], Fcart(th + h)[1])
    return ((Pp[0] - Fm[0]) / (2 * h), (Pp[1] - Fm[1]) / (2 * h))

# compatibility check of displacement construction: V_th' + V_r = g_thth
thck = np.linspace(0.01, np.pi - 0.01, 999)
Vr, Vt, grr, dgrr, grt = Vfun(thck)
dVt = np.gradient(Vt, thck)
g = ang(thck)
print(f"compatibility max|V_th'+V_r-g_tt| = {np.max(np.abs(dVt + Vr - g['gtt'])):.3e}")

# strain check: strains from u vs tangent law at several angles
print("\nstrain check (u from closed form vs tangent law), r=0.003:")
r0 = 3e-3
for thd in [0.0, 30.0, 45.0, 60.0, 90.0, 120.0]:
    thv = np.deg2rad(thd)
    F = Fcart(thv); Fp = Fcart_d(thv)
    exx = phi0 / np.sqrt(r0) * (0.5 * np.cos(thv) * F[0] - np.sin(thv) * Fp[0])
    eyy = phi0 / np.sqrt(r0) * (0.5 * np.sin(thv) * F[1] + np.cos(thv) * Fp[1])
    exy = 0.5 * phi0 / np.sqrt(r0) * (0.5 * np.sin(thv) * F[0] + np.cos(thv) * Fp[0]
                                      + 0.5 * np.cos(thv) * F[1] - np.sin(thv) * Fp[1])
    gg = ang(thv)
    lx, ly, lxy = phi0 / np.sqrt(r0) * gg["gxx"], phi0 / np.sqrt(r0) * gg["gyy"], phi0 / np.sqrt(r0) * gg["gxy"]
    print(f"  th={thd:5.1f}: e_xx {exx:+.6e} vs {lx:+.6e} | e_yy {eyy:+.6e} vs {ly:+.6e} | e_xy {exy:+.6e} vs {lxy:+.6e}")

# C_y and CTOD checks
F90 = Fcart(np.pi / 2)
print(f"\nu_y(r,90)/sqrt(r) = phi0*F_y(90) = {phi0*F90[1]:.8e} (target C_y)")
print(f"2*C_yy(90)         = {2*phi0*float(ang(np.pi/2)['gyy']):.8e}")
uxp, uyp = Fcart(np.pi)
print(f"crack faces theta=+-pi: u_x = {phi0*uxp:.4e} (0), u_y = {phi0*uyp:.4e} >0; CTOD at r=1: {2*phi0*uyp:.6e} vs textbook 8K/E*sqrt(1/2pi) = {8*K/Et/np.sqrt(2*np.pi):.6e}")

# ---------------- J contour ----------------
def J_contour(rho, mode):
    thg = np.linspace(-np.pi + 1e-9, np.pi - 1e-9, 400001)
    q = K / np.sqrt(2 * np.pi * rho)
    gg = ang(thg)
    sxx, syy, sxy = q * gg["fxx"], q * gg["fyy"], q * gg["fxy"]
    F = Fcart(thg); Fp = Fcart_d(thg)
    ct, st = np.cos(thg), np.sin(thg)
    uxpx = phi0 / np.sqrt(rho) * (0.5 * ct * F[0] - st * Fp[0])
    uypx = phi0 / np.sqrt(rho) * (0.5 * ct * F[1] - st * Fp[1])
    tx, ty = sxx * ct + sxy * st, sxy * ct + syy * st
    if mode == "lin":
        exx = (sxx - nt * syy) / Et; eyy = (syy - nt * sxx) / Et; exy = (1 + nt) * sxy / Et
        W = 0.5 * (sxx * exx + syy * eyy + 2 * sxy * exy)
    else:
        sdot = sxx**2 + syy**2 + 2 * sxy**2
        skk = sxx + syy
        ss = sdot - skk**2 / 3.0
        se = np.sqrt(1.5 * ss)
        tym = np.minimum(np.where(se > 0, sy / np.maximum(se, 1e-300), 1.0), 1.0)
        W = 0.5 * ((1 + nu) / E * sdot - nu / E * skk**2) + 0.75 * (1 / Et - 1 / E) * ss * (1 - tym**2)
    return trapz(W * ct - (tx * uxpx + ty * uypx), thg) * rho

print("\nJ contour (expect 0.800 N/mm):")
for rho in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    jl = J_contour(rho, "lin"); jf = J_contour(rho, "full")
    print(f"  rho={rho:.0e}: J(tangent W)={jl:.8f}, J(full-law W)={jf:.6f} (err {jf-J:+.2e})")
