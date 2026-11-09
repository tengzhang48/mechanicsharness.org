import numpy as np
from scipy.optimize import brentq

E = 70000.0; nu = 0.30; sy = 10.0; Jval = 0.800; Czz_meas = -0.005

# ---- effective moduli of the asymptotic (sigma_e -> inf) reduction ----
# eps = (1+nu)/E sig - nu/E sig_kk delta + (3/2)A s,  A = 1/T - 1/E
# == (1+nu*)/E* sig - nu*/E* sig_kk delta  with E* = T, nu* = nu*t + (1-t)/2, t = T/E
def nup_of_t(t): return nu*t + (1.0 - t)/2.0

def Czz_pred(t):
    T = t*E; nup = nup_of_t(t); K = np.sqrt(Jval*T)
    return -2.0*nup*K/(T*np.sqrt(2.0*np.pi))   # eps_zz(r,0) = -(2 nu*/T) K/sqrt(2 pi r)

f = lambda t: Czz_pred(t) - Czz_meas
t = brentq(f, 1e-12, 1.0-1e-12, xtol=1e-15, rtol=8.9e-16, maxiter=200)
T = t*E; nup = nup_of_t(t); K = np.sqrt(Jval*T)
print(f"t = T/E = {t:.8f}")
print(f"1) T = {T:.2f} MPa   nu* = {nup:.8f}   K = {K:.5f} MPa sqrt(mm)")
print(f"   C_zz reproduced = {Czz_pred(t):.8f} (measured {Czz_meas})")

# ---- leading field: plane-stress Williams mode I with (T, nu*) ----
def stress(r, th):
    s = K/np.sqrt(2*np.pi*r)
    c2, s2 = np.cos(th/2), np.sin(th/2)
    s3, c3 = np.sin(3*th/2), np.cos(3*th/2)
    return s*c2*(1-s2*s3), s*c2*(1+s2*s3), s*c2*s2*c3  # sxx, syy, sxy

def strain_cart(r, th):
    sxx, syy, sxy = stress(r, th)
    return (sxx-nup*syy)/T, (syy-nup*sxx)/T, (1+nup)*sxy/T

# ---- part 2: theta = 60 deg ----
th = np.deg2rad(60.0); c, s = np.cos(th), np.sin(th)
sxx, syy, sxy = stress(1.0, th)
exx, eyy, exy = strain_cart(1.0, th)
srr = sxx*c*c + syy*s*s + 2*sxy*s*c
stt = sxx*s*s + syy*c*c - 2*sxy*s*c
srt = (syy-sxx)*s*c + sxy*(c*c - s*s)
err = exx*c*c + eyy*s*s + 2*exy*s*c
ett = exx*s*s + eyy*c*c - 2*exy*s*c
ert = (eyy-exx)*s*c + exy*(c*c - s*s)
print("\n2) theta=60deg:  (polar-transform of cartesian strains)")
print(f"   C_rr    = {err:.6e} mm^1/2   (direct (srr-nu*stt)/T: {(srr-nup*stt)/T:.6e})")
print(f"   C_thth  = {ett:.6e} mm^1/2   (direct:               {(stt-nup*srr)/T:.6e})")
print(f"   C_grth  = {2*ert:.6e} mm^1/2 (sigma_rth at 60deg = {srt:.3e})")

# ---- part 3: u_y(r,90deg)/sqrt(r) ----
kappa = (3-nup)/(1+nup)
def disp(r, th):
    f = K*(1+nup)/T*np.sqrt(r/(2*np.pi))
    c2, s2 = np.cos(th/2), np.sin(th/2)
    return f*c2*(kappa-1+2*s2*s2), f*s2*(kappa+1-2*c2*c2)
uy = disp(1.0, np.pi/2)[1]
print(f"\n3) C_y(90deg) = {uy:.6e} mm^1/2   (kappa* = {kappa:.6f})")

# ---- part 4: screening interval on theta = 0 (leading field) ----
# sigma_e(r,0) = K/sqrt(2 pi r) (plane stress, equibiaxial); principal strains:
# eps_xx=eps_yy=(1-nu*)sigma/T, eps_zz=-2 nu* sigma/T -> max |.| = |eps_zz|
sig = lambda r: K/np.sqrt(2*np.pi*r)
r_hi_se = (K/(10*sy))**2/(2*np.pi)
epszz = lambda r: -2*nup*sig(r)/T
r_lo_st = brentq(lambda r: abs(epszz(r)) - 0.05, 1e-12, 1.0)
print(f"\n4) screens on theta=0:")
print(f"   sigma_e >= 100 MPa  -> r <= {r_hi_se:.6f} mm")
print(f"   |eps_zz| <= 0.05    -> r >= {r_lo_st:.6f} mm  (eps_zz(r,0)=-0.005/sqrt(r): {(abs(epszz(r_lo_st))*np.sqrt(r_lo_st)):.6f})")
print(f"   r <= 0.200 mm       -> r <= 0.200000 mm")
r_lo, r_hi = r_lo_st, min(r_hi_se, 0.200)
print(f"   interval = [{r_lo:.6f}, {r_hi:.6f}] mm, empty? {r_hi < r_lo}")

# ---- verification A: numerical J integral on the leading field ----
def u_cart(x, y):
    r = np.hypot(x, y); th = np.arctan2(y, x)
    return disp(r, th)
def J_integral(R, n=400001):
    th = np.linspace(-np.pi, np.pi, n)
    r = R
    sxx, syy, sxy = stress(r, th)
    exx, eyy, exy = strain_cart(r, th)
    W = 0.5*(sxx*exx + syy*eyy + 2*sxy*exy)
    cth, sth = np.cos(th), np.sin(th)
    tx, ty = sxx*cth + sxy*sth, sxy*cth + syy*sth
    x, y = R*cth, R*sth; h = 1e-6*R
    ux_x = (u_cart(x+h, y)[0] - u_cart(x-h, y)[0])/(2*h)
    uy_x = (u_cart(x+h, y)[1] - u_cart(x-h, y)[1])/(2*h)
    integrand = W*cth - (tx*ux_x + ty*uy_x)
    return np.trapezoid(integrand, th)*R
print("\nverif A: J from contour integral on leading field vs K^2/T:")
for R in (1e-3, 1e-2, 1e-1):
    print(f"   R={R:7.4f} mm -> J = {J_integral(R):.8f}  (K^2/T = {K*K/T:.8f}, given J = {Jval})")

# ---- verification B: strains from finite differences of u ----
print("\nverif B: FD strains from displacement field vs stress-derived:")
h = 1e-6
for (r0, th0) in [(0.005, 0.6), (0.02, 1.1), (0.05, 2.0), (0.1, -1.3)]:
    x0, y0 = r0*np.cos(th0), r0*np.sin(th0)
    exx_fd = (u_cart(x0+h, y0)[0]-u_cart(x0-h, y0)[0])/(2*h)
    eyy_fd = (u_cart(x0, y0+h)[1]-u_cart(x0, y0-h)[1])/(2*h)
    exy_fd = 0.5*((u_cart(x0, y0+h)[0]-u_cart(x0, y0-h)[0])/(2*h)
                + (u_cart(x0+h, y0)[1]-u_cart(x0-h, y0)[1])/(2*h))
    ean = strain_cart(r0, th0)
    print(f"   r={r0:6.3f} th={th0:+5.2f}: FD {exx_fd:+.6e} {eyy_fd:+.6e} {exy_fd:+.6e}"
          f" | analytic {ean[0]:+.6e} {ean[1]:+.6e} {ean[2]:+.6e}")

# ---- verification C: consistency of the true law on theta=0 in the annulus ----
print("\nverif C: true vs leading eps_zz on theta=0 (plastic-reduction error):")
A = 1/T - 1/E
for r0 in (r_lo, (r_lo+r_hi)/2, r_hi):
    sg = sig(r0); se = sg
    e_true = -2*sg*(nu/E + 0.5*A*(1 - sy/se))
    e_lead = -2*sg*(nu/E + 0.5*A)
    print(f"   r={r0:.6f}: sigma_e={se:8.3f} MPa, rel.diff = {abs((e_true-e_lead)/e_lead)*100:.2f}% (sigma_y/sigma_e={sy/se:.3f})")
