import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

trapz = getattr(np, "trapezoid", None) or np.trapz

E, nu, sy, J = 70000.0, 0.30, 10.0, 0.800
Czz_meas = -0.005
Tsol = brentq(lambda T: -(0.5 - (0.5 - nu) * T / E) / T * 2.0 * np.sqrt(J * T) / np.sqrt(2 * np.pi) - Czz_meas,
              100.0, 69999.0, xtol=1e-13, rtol=1e-15)
Et, nt, K = Tsol, 0.5 - (0.5 - nu) * Tsol / E, np.sqrt(J * Tsol)
phi0 = K / (Et * np.sqrt(2 * np.pi))
print(f"T={Tsol:.6f} MPa, E_t={Et:.4f}, nu_t={nt:.8f}, K_t={K:.6f}, phi0={phi0:.8e}")

def ang(th):
    th = np.atleast_1d(np.asarray(th, float))
    c, s = np.cos(th / 2), np.sin(th / 2)
    s3, c3 = np.sin(1.5 * th), np.cos(1.5 * th)
    fxx = c * (1 - s * s3); fyy = c * (1 + s * s3); fxy = c * s * c3
    dfxx = -(s / 2) * (1 - s * s3) + c * (-(c / 2) * s3 - s * 1.5 * c3)
    dfyy = -(s / 2) * (1 + s * s3) + c * ((c / 2) * s3 + s * 1.5 * c3)
    dfxy = 0.5 * (c * c - s * s) * c3 - 1.5 * c * s * s3
    gxx, gyy, gxy = fxx - nt * fyy, fyy - nt * fxx, (1 + nt) * fxy
    dgxx, dgyy, dgxy = dfxx - nt * dfyy, dfyy - nt * dfxx, (1 + nt) * dfxy
    srr, stt, srt = c * (1 + s * s), c * (1 - s * s), s * c * c
    return dict(fxx=fxx, fyy=fyy, fxy=fxy, gxx=gxx, gyy=gyy, gxy=gxy,
                dgxx=dgxx, dgyy=dgyy, dgxy=dgxy, srr=srr, stt=stt, srt=srt,
                gtt=stt - nt * srr, grt=(1 + nt) * srt, th=th)

def rhs(th, y):
    g = ang(th)
    st, ct = np.sin(th), np.cos(th)
    gxx, gyy, gxy = g["gxx"][0], g["gyy"][0], g["gxy"][0]
    dgxx, dgyy, dgxy = g["dgxx"][0], g["dgyy"][0], g["dgxy"][0]
    Awx = -0.5 * st * gxx + ct * dgxx + 0.5 * ct * gxy + st * dgxy
    Awy = -0.5 * st * gxy + ct * dgxy + 0.5 * ct * gyy + st * dgyy
    w = y[0]
    wt = -st * Awx + ct * Awy          # d(w~)/d(theta)
    grt, gtt = g["grt"][0], g["gtt"][0]
    er, et = np.array([ct, st]), np.array([-st, ct])
    Fp = (grt - w) * er + gtt * et     # dF/d(theta)
    return [wt, Fp[0], Fp[1]]

sol = solve_ivp(rhs, (0.0, np.pi), [0.0, 2.0 * (1.0 - nt), 0.0],
                dense_output=True, rtol=1e-12, atol=1e-14, max_step=1e-3)
ths = np.linspace(0.0, np.pi, 400001)
Y = sol.sol(ths)
wt_, Fx, Fy = Y[0], Y[1], Y[2]
spFx, spFy = CubicSpline(ths, Fx), CubicSpline(ths, Fy)
spFx1, spFy1 = spFx.derivative(), spFy.derivative()

# consistency check of the two rotation-gradient equations
gA = ang(ths)
Awx = -0.5 * np.sin(ths) * gA["gxx"] + np.cos(ths) * gA["dgxx"] + 0.5 * np.cos(ths) * gA["gxy"] + np.sin(ths) * gA["dgxy"]
Awy = -0.5 * np.sin(ths) * gA["gxy"] + np.cos(ths) * gA["dgxy"] + 0.5 * np.cos(ths) * gA["gyy"] + np.sin(ths) * gA["dgyy"]
resid_A = np.cos(ths) * Awx + np.sin(ths) * Awy + 0.5 * wt_
print(f"rotation-ODE cross-check max|residual| = {np.max(np.abs(resid_A)):.3e} (should be ~0)")

# F_y(pi/2) vs required (3-nu_t)/sqrt(2) = 2*g_yy(90)
print(f"F_y(pi/2) = {spFy(np.pi/2):.10f}, expected (3-nu_t)/sqrt2 = {(3-nt)/np.sqrt(2):.10f}")
print(f"C_y = phi0*F_y(pi/2) = {phi0*spFy(np.pi/2):.8e} mm^0.5")
print(f"2*C_yy(90) = {2*phi0*(np.cos(np.pi/4)*(1+np.sin(np.pi/4)*np.sin(3*np.pi/4)) - nt*np.cos(np.pi/4)*(1-np.sin(np.pi/4)*np.sin(3*np.pi/4))):.8e}")

# strain check from reconstructed u at various angles (incl. shear != 0)
def u_and_grad(r, th):
    a = np.abs(th); sgn = np.sign(np.th)if False else np.sign(th)
    Fxv = spFx(a); Fyv = spFy(a) * sgn
    Fxp = spFx1(a); Fyp = spFy1(a) * sgn
    F = np.array([Fxv, Fyv]); Fp = np.array([Fxp, Fyp])
    u = phi0 * np.sqrt(r) * F
    ux_x = phi0 / np.sqrt(r) * (0.5 * np.cos(th) * F[0] - np.sin(th) * Fp[0])
    uy_x = phi0 / np.sqrt(r) * (0.5 * np.cos(th) * F[1] - np.sin(th) * Fp[1])
    ux_y = phi0 / np.sqrt(r) * (0.5 * np.sin(th) * F[0] + np.cos(th) * Fp[0])
    uy_y = phi0 / np.sqrt(r) * (0.5 * np.sin(th) * F[1] + np.cos(th) * Fp[1])
    return u, ux_x, ux_y, uy_x, uy_y

print("\nstrain check (reconstructed u vs tangent law), r=0.003:")
r0 = 3e-3
for thd in [0.0, 30.0, 45.0, 60.0, 90.0]:
    thv = np.deg2rad(thd)
    u, uxx, uxy, uyx, uyy = u_and_grad(r0, thv)
    g = ang(thv)
    law_xx, law_yy, law_xy = phi0 / np.sqrt(r0) * g["gxx"][0], phi0 / np.sqrt(r0) * g["gyy"][0], phi0 / np.sqrt(r0) * g["gxy"][0]
    print(f"  th={thd:5.1f}: e_xx {uxx:.8e} vs {law_xx:.8e} | e_yy {uyy:.8e} vs {law_yy:.8e} | e_xy {0.5*(uxy+uyx):+.8e} vs {law_xy:+.8e}")

# compare with the (wrong) textbook-recalled formula
def old_G(th):
    c, s = np.cos(th / 2), np.sin(th / 2)
    return 2 * c * ((1 - nt) + (1 + nt) * s**2), 2 * s * (2 - (1 + nt) * c**2)
for thd in [30.0, 60.0, 90.0]:
    thv = np.deg2rad(thd)
    Gx, Gy = old_G(thv)
    print(f"  th={thd:5.1f}: old formula G=({Gx:+.6f},{Gy:+.6f})  reconstructed F=({spFx(thv):+.6f},{spFy(thv):+.6f})")

# ---------------- J contour with reconstructed u ----------------
def J_contour(rho, mode):
    thg = np.linspace(-np.pi + 1e-12, np.pi - 1e-12, 400001)
    q = K / np.sqrt(2 * np.pi * rho)
    g = ang(thg)
    sxx, syy, sxy = q * g["fxx"], q * g["fyy"], q * g["fxy"]
    sgn = np.sign(thg)
    F = np.array([spFx(np.abs(thg)), spFy(np.abs(thg)) * sgn])
    Fp = np.array([spFx1(np.abs(thg)), spFy1(np.abs(thg)) * sgn])
    ct, st = np.cos(thg), np.sin(thg)
    ux = phi0 * np.sqrt(rho) * F[0]; uy = phi0 * np.sqrt(rho) * F[1]
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
