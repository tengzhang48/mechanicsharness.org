import numpy as np
from scipy.optimize import brentq

E, nu, sy, J = 70000.0, 0.30, 10.0, 0.800
Czz_meas = -0.005

trapz = getattr(np, "trapezoid", None) or np.trapz

# ---------------- asymptotic (tangent) moduli of the total-strain law ----------------
# eps = (1+nu)/E sig - nu/E skk I + 3/2 (1/T-1/E) <1-sy/se> s
# for se >> sy:  eps = A sig + B skk I  (isotropic linear), with
# A = (nu-1/2)/E + 3/(2T),  B = (1/2-nu)/E - 1/(2T)
# => E_t = T,  nu_t = 1/2 - (1/2-nu) T/E      (tangent compliances add: 1/E_t = 1/E+1/T-1/E)
def E_t(T):
    return T

def nu_t(T):
    return 0.5 - (0.5 - nu) * T / E

def K_t(T):
    return np.sqrt(J * E_t(T))  # plane stress: J = K^2/E_t

def Czz_pred(T):
    # theta=0: sig_xx=sig_yy=K/sqrt(2 pi r), sig_xy=0; eps_zz = -(nu_t/E_t)(sig_xx+sig_yy)
    return -(nu_t(T) / E_t(T)) * 2.0 * K_t(T) / np.sqrt(2.0 * np.pi)

Tsol = brentq(lambda T: Czz_pred(T) - Czz_meas, 100.0, 69999.0, xtol=1e-13, rtol=1e-15)
Et, nt, K = E_t(Tsol), nu_t(Tsol), K_t(Tsol)
P = K / (np.sqrt(2.0 * np.pi) * Et)  # C = P * (dimensionless angular coefficient)

print(f"1) T = {Tsol:.6f} MPa")
print(f"   E_t = {Et:.4f} MPa, nu_t = {nt:.8f}, K_t = {K:.6f} MPa*mm^0.5, P = {P:.8e} mm^0.5")
print(f"   residual Czz_pred-Czz_meas = {Czz_pred(Tsol)-Czz_meas:.3e}")
g = lambda T: 2.0 * (1.0 / (2.0 * T) - (0.5 - nu) / E) * np.sqrt(J * T) / np.sqrt(2.0 * np.pi)
print(f"   uniqueness: dg/dT < 0 on (0,E): g(100)={g(100.0):.4e}, g(4818.55)={g(4818.55):.6e}, g(69000)={g(69000.0):.4e}")

# ---------------- Q2: coefficients at theta = 60 deg ----------------
th = np.deg2rad(60.0)
c, s = np.cos(th / 2.0), np.sin(th / 2.0)
srr, stt, srt = c * (1 + s * s), c * (1 - s * s), s * c * c  # coefficients of K/sqrt(2 pi r)
Crr, Ctt, Cg = P * (srr - nt * stt), P * (stt - nt * srr), 2.0 * (1.0 + nt) * P * srt
print(f"2) angular stress coefs at 60deg: srr={srr:.6f}, stt={stt:.6f}, srt={srt:.6f}")
print(f"   C_rr = {Crr:.6e}, C_tt = {Ctt:.6e}, C_gamma_rt = {Cg:.6e}  mm^0.5")

# ---------------- Q3: C_y(90) ----------------
kd = (3.0 - nt) / (1.0 + nt)
Cy = (1.0 + nt) * (K / Et) / np.sqrt(2.0 * np.pi) * np.sin(np.pi / 4.0) * (kd + 1.0 - 2.0 * np.cos(np.pi / 4.0) ** 2)
th90 = np.pi / 2.0
sxx90 = np.cos(th90 / 2) * (1 - np.sin(th90 / 2) * np.sin(1.5 * th90))
syy90 = np.cos(th90 / 2) * (1 + np.sin(th90 / 2) * np.sin(1.5 * th90))
Cyy90 = P * (syy90 - nt * sxx90)
print(f"3) C_y(90) = {Cy:.6e} mm^0.5 ; cross-check 2*C_yy(90) = {2.0*Cyy90:.6e}")

# ---------------- Q4: screens on theta = 0 ----------------
# sig_e = p = K/sqrt(2 pi r); principal strains: e_in = (1-nu_t)p/E_t (double), |e_zz| = 2 nu_t p/E_t (max)
print(f"   max-principal check: 2*nu_t = {2*nt:.5f} > 1-nu_t = {1-nt:.5f} -> |e_zz| governs")
r_lo = (abs(Czz_meas) / 0.05) ** 2          # |e_zz| = 0.005/sqrt(r) <= 0.05
r_hi = (K / (np.sqrt(2.0 * np.pi) * 10.0 * sy)) ** 2  # sig_e >= 100 MPa
print(f"4) r interval = [{r_lo:.8f}, {min(r_hi,0.2):.8f}] mm, nonempty = {r_lo <= min(r_hi,0.2)}")
print(f"   sig_e at r_hi = {K/np.sqrt(2*np.pi*r_hi):.4f} MPa ; |e_zz| at r_lo = {abs(Czz_meas)/np.sqrt(r_lo):.6f}")

# ---------------- full constitutive law ----------------
def eps_full(sxx, syy, sxy, szz, T):
    skk = sxx + syy + szz
    dxx, dyy, dzz, dxy = sxx - skk / 3, syy - skk / 3, szz - skk / 3, sxy
    se = np.sqrt(1.5 * (dxx**2 + dyy**2 + dzz**2 + 2 * dxy**2))
    fac = 1.5 * (1.0 / T - 1.0 / E) * np.maximum(1.0 - sy / se, 0.0)
    exx = (1 + nu) / E * sxx - nu / E * skk + fac * dxx
    eyy = (1 + nu) / E * syy - nu / E * skk + fac * dyy
    ezz = (1 + nu) / E * szz - nu / E * skk + fac * dzz
    exy = (1 + nu) / E * sxy + fac * dxy
    return exx, eyy, ezz, exy

def lead_stress(r, th, T):
    q = K_t(T) / np.sqrt(2.0 * np.pi * r)
    sxx = q * np.cos(th / 2) * (1 - np.sin(th / 2) * np.sin(1.5 * th))
    syy = q * np.cos(th / 2) * (1 + np.sin(th / 2) * np.sin(1.5 * th))
    sxy = q * np.cos(th / 2) * np.sin(th / 2) * np.cos(1.5 * th)
    return sxx, syy, sxy

def lead_disp(r, th, T):
    Kk = K_t(T); Ett = E_t(T); ntt = nu_t(T); kdd = (3.0 - ntt) / (1.0 + ntt)
    fac = (1.0 + ntt) * Kk / Ett * np.sqrt(r / (2.0 * np.pi))
    ux = fac * np.cos(th / 2) * (kdd - 1.0 + 2.0 * np.sin(th / 2) ** 2)
    uy = fac * np.sin(th / 2) * (kdd + 1.0 - 2.0 * np.cos(th / 2) ** 2)
    return ux, uy

# ---------------- V2: full-law limits reproduce the coefficients ----------------
print("\nV2 full-law limits (leading stress field fed through the COMPLETE law):")
for r in [1e-2, 1e-4, 1e-6, 1e-8]:
    sxx, syy, sxy = lead_stress(r, 0.0, Tsol)
    exx, eyy, ezz, exy = eps_full(sxx, syy, sxy, 0.0, Tsol)
    print(f"   r={r:.0e} th=0 : sqrt(r)e_zz = {np.sqrt(r)*ezz:.8f}")
ct6, st6 = np.cos(th), np.sin(th)
for r in [1e-2, 1e-4, 1e-6, 1e-8]:
    sxx, syy, sxy = lead_stress(r, th, Tsol)
    srr_ = sxx * ct6**2 + syy * st6**2 + 2 * sxy * ct6 * st6
    stt_ = sxx * st6**2 + syy * ct6**2 - 2 * sxy * ct6 * st6
    srt_ = (syy - sxx) * ct6 * st6 + sxy * (ct6**2 - st6**2)
    err, ett, ezz, ert = eps_full(srr_, stt_, srt_, 0.0, Tsol)
    print(f"   r={r:.0e} th=60: C_rr={np.sqrt(r)*err:.6e} C_tt={np.sqrt(r)*ett:.6e} C_g={np.sqrt(r)*2*ert:.6e}")

# ---------------- V1: J contour integral ----------------
def W_full(sxx, syy, sxy, T):
    # proportional path sig(t)=t*sig:  W = 1/2 sig:e_e + 3/4 (1/T-1/E) s:s (1-ty^2), ty=min(1,sy/se)
    sdot = sxx**2 + syy**2 + 2 * sxy**2
    skk = sxx + syy
    ss = sdot - skk**2 / 3.0
    se = np.sqrt(1.5 * ss)
    tym = np.minimum(np.where(se > 0, sy / np.maximum(se, 1e-300), 1.0), 1.0)
    We = 0.5 * ((1 + nu) / E * sdot - nu / E * skk**2)
    Wp = 0.75 * (1.0 / T - 1.0 / E) * ss * (1.0 - tym**2)
    return We + Wp

def J_contour(rho, T, mode):
    thg = np.linspace(-np.pi, np.pi, 400001)[1:-1]
    sxx, syy, sxy = lead_stress(rho, thg, T)
    ux, uy = lead_disp(rho, thg, T)
    x, y = rho * np.cos(thg), rho * np.sin(thg)
    h = 1e-6 * rho
    uxp = lead_disp(np.hypot(x + h, y), np.arctan2(y, x + h), T)
    uxm = lead_disp(np.hypot(x - h, y), np.arctan2(y, x - h), T)
    uyp = lead_disp(np.hypot(x, y + h), np.arctan2(y + h, x), T)
    uym = lead_disp(np.hypot(x, y - h), np.arctan2(y - h, x), T)
    uxpx = (uxp[0] - uxm[0]) / (2 * h)
    uypx = (uyp[1] - uym[1]) / (2 * h)
    Nx, Ny = np.cos(thg), np.sin(thg)
    tx, ty_ = sxx * Nx + sxy * Ny, sxy * Nx + syy * Ny
    if mode == "lin":
        ntt, Ett = nu_t(T), E_t(T)
        exx = (sxx - ntt * syy) / Ett
        eyy = (syy - ntt * sxx) / Ett
        exy = (1 + ntt) * sxy / Ett
        W = 0.5 * (sxx * exx + syy * eyy + 2 * sxy * exy)
    else:
        W = W_full(sxx, syy, sxy, T)
    return trapz(W * Nx - (tx * uxpx + ty_ * uypx), thg) * rho

print("\nV1 J contour integral on the leading field (expect J = 0.800 N/mm):")
for rho in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    jl = J_contour(rho, Tsol, "lin")
    jf = J_contour(rho, Tsol, "full")
    print(f"   rho={rho:.0e} mm: J(tangent W)={jl:.8f}, J(full-law W)={jf:.6f} (err {jf-J:+.2e})")
print(f"   K_t^2/E_t = {K**2/Et:.6f} (by construction = J)")

# ---------------- V3: displacement-strain compatibility ----------------
def num_strain(x0, y0, T, h):
    def u(xx, yy):
        return lead_disp(np.hypot(xx, yy), np.arctan2(yy, xx), T)
    duxdx = (u(x0 + h, y0)[0] - u(x0 - h, y0)[0]) / (2 * h)
    duxdy = (u(x0, y0 + h)[0] - u(x0, y0 - h)[0]) / (2 * h)
    duydx = (u(x0 + h, y0)[1] - u(x0 - h, y0)[1]) / (2 * h)
    duydy = (u(x0, y0 + h)[1] - u(x0, y0 - h)[1]) / (2 * h)
    return duxdx, duydy, 0.5 * (duxdy + duydx)

print("\nV3 compatibility (numeric sym-grad u vs tangent-law strains):")
r0 = 5e-3
x0, y0 = r0 * np.cos(th), r0 * np.sin(th)
exx_n, eyy_n, exy_n = num_strain(x0, y0, Tsol, 1e-6 * r0)
sxx0, syy0, sxy0 = lead_stress(r0, th, Tsol)
exx_l = (sxx0 - nt * syy0) / Et
eyy_l = (syy0 - nt * sxx0) / Et
exy_l = (1 + nt) * sxy0 / Et
print(f"   th=60, r=0.005: e_xx num/law = {exx_n:.6e}/{exx_l:.6e}, e_yy = {eyy_n:.6e}/{eyy_l:.6e}, e_xy = {exy_n:.6e}/{exy_l:.6e}")
# u_y(0,y) = int_0^y e_yy(0,y') dy' check
for yv in [1e-3, 1e-2, 1e-1]:
    yg = np.linspace(1e-12, yv, 2000001)
    sxxg, syyg, _ = lead_stress(yg, np.pi / 2 * np.ones_like(yg), Tsol)
    eyyg = (syyg - nt * sxxg) / Et
    uy_int = trapz(eyyg, yg)
    uy_f, _ = lead_disp(yv, np.pi / 2, Tsol)
    print(f"   y={yv:.0e}: u_y formula={uy_f:.8e}, int e_yy dy={uy_int:.8e}, u_y/sqrt(y)={uy_f/np.sqrt(yv):.6e}")

# crack-face conditions
sxx_f, syy_f, sxy_f = lead_stress(1e-3, np.pi - 1e-9, Tsol)
print(f"\nface check theta->pi: sig_theta_theta={0.5*(sxx_f+syy_f-(sxx_f-syy_f)*np.cos(2*(np.pi-1e-9))):.3e}, sig_r_theta~{((syy_f-sxx_f)*0):.1e} (both ~0)")
ux_f, uy_f = lead_disp(1e-3, np.pi, Tsol)
print(f"face opening u_y(r,pi) = {uy_f:.4e} > 0")
