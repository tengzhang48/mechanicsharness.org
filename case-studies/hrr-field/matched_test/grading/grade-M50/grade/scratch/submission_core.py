"""Infer tangent modulus T from C_zz(0), then compute requested coefficients.

Model: total-strain (deformation-theory) bilinear law. As r->0, sigma_e -> infinity, so
<1 - sigma_y/sigma_e> -> 1 and the law becomes, to leading order, a linear isotropic law:
    eps_ij = [(1+nu)/E + 1.5(1/T - 1/E)] sig_ij + [-nu/E - 0.5(1/T - 1/E)] sig_kk delta_ij
i.e. Hooke with E* = T and nu* = 1/2 + (nu - 1/2) T/E (plane stress, sigma_zz = 0).
Leading field: mode I Williams field, p = 1/2, K = sqrt(J E*) = sqrt(J T).
Units: mm, N, MPa.
"""
import numpy as np
from scipy.optimize import brentq

E, nu, sy, J = 70000.0, 0.30, 10.0, 0.800
CZZ_MEAS = -0.005

nustar = lambda T: 0.5 + (nu - 0.5) * T / E
Kof = lambda T: np.sqrt(J * T)

def Czz_lead(T):
    return -2.0 * nustar(T) * np.sqrt(J / T) / np.sqrt(2.0 * np.pi)

# ---- 1. solve for T --------------------------------------------------------
T = brentq(lambda t: Czz_lead(t) - CZZ_MEAS, 100.0, E - 100.0, xtol=1e-13, rtol=8.9e-16)
nus, K = nustar(T), Kof(T)
A = K / np.sqrt(2.0 * np.pi)          # angular amplitude: sig_ij = A r^-1/2 S_ij(th)
print(f"1) T  = {T:.6f} MPa   (nu* = {nus:.6f}, K = {K:.6f} MPa mm^1/2)")
print(f"   residual |Czz_lead + 0.005| = {abs(Czz_lead(T) - CZZ_MEAS):.2e}")

# ---- full constitutive law (for verification) ------------------------------
def full_strain(sxx, syy, sxy, T):
    """Plane stress (szz = 0). Returns eps_xx, eps_yy, eps_xy(tensor), eps_zz, sigma_e."""
    sig = np.array([sxx, syy, 0.0, sxy])
    m = (sxx + syy) / 3.0
    s = np.array([sxx - m, syy - m, -m, sxy])
    se = np.sqrt(1.5 * (s[0] ** 2 + s[1] ** 2 + s[3] ** 2 + 2 * s[2] ** 2 + 2 * s[2] ** 2 * 0
                        + s[2] ** 2 * 0 + s[0] * 0) )
    # careful: components (xx,yy,zz,xy); s_ij s_ij sums diagonals + 2*offdiag
    se = np.sqrt(1.5 * (s[0] ** 2 + s[1] ** 2 + s[2] ** 2 + 2 * s[3] ** 2))
    f = 1.5 * (1.0 / T - 1.0 / E) * np.maximum(1.0 - sy / se, 0.0)
    eel = np.array([(1 + nu) / E * sig[i] for i in range(4)])
    eel[0] -= nu / E * (sxx + syy); eel[1] -= nu / E * (sxx + syy); eel[2] -= nu / E * (sxx + syy)
    eps = eel + f * s
    return eps[0], eps[1], eps[3], eps[2], se

# K-field stress (Cartesian), mode I
def sig_K(r, th, K):
    f = K / np.sqrt(2 * np.pi * r)
    c2, s2 = np.cos(th / 2), np.sin(th / 2)
    sxx = f * c2 * (1 - s2 * np.sin(3 * th / 2))
    syy = f * c2 * (1 + s2 * np.sin(3 * th / 2))
    sxy = f * c2 * s2 * np.cos(3 * th / 2)
    return sxx, syy, sxy

def cart2pol_strain(exx, eyy, exy, th):
    c, s = np.cos(th), np.sin(th)
    err = exx * c**2 + eyy * s**2 + 2 * exy * s * c
    ett = exx * s**2 + eyy * c**2 - 2 * exy * s * c
    ert = (eyy - exx) * s * c + exy * (c**2 - s**2)
    return err, ett, ert

# ---- 2. coefficients at theta = 60 deg (full law, r -> 0 extrapolation) ----
th = np.deg2rad(60.0)
print("\n2) theta = 60 deg, sqrt(r)*coefficient from FULL law vs leading field:")
print("   r        C_rr        C_tt       C_grt      (dev of e_zz-like check omitted)")
for r in (0.2, 0.05, 0.01, 0.002, 0.0005):
    sxx, syy, sxy = sig_K(r, th, K)
    exx, eyy, exy, ezz, se = full_strain(sxx, syy, sxy, T)
    err, ett, ert = cart2pol_strain(exx, eyy, exy, th)
    print(f"   {r:7.4f}  {np.sqrt(r)*err: .6f}  {np.sqrt(r)*ett: .6f}  {np.sqrt(r)*2*ert: .6f}")

# leading values (effective linear, closed form)
def lead_coeffs(th, T, K):
    nus = nustar(T)
    # polar stresses at any theta from Cartesian K-field, then plane-stress Hooke with (T,nu*)
    def S(th):
        c2, s2 = np.cos(th/2), np.sin(th/2)
        Sxx = c2 * (1 - s2*np.sin(3*th/2)); Syy = c2 * (1 + s2*np.sin(3*th/2)); Sxy = c2*s2*np.cos(3*th/2)
        c, s = np.cos(th), np.sin(th)
        Srr = Sxx*c**2 + Syy*s**2 + 2*Sxy*s*c
        Stt = Sxx*s**2 + Syy*c**2 - 2*Sxy*s*c
        Srt = (Syy - Sxx)*s*c + Sxy*(c**2 - s**2)
        return Srr, Stt, Srt
    Srr, Stt, Srt = S(th)
    a = K / (T * np.sqrt(2*np.pi))
    return a*(Srr - nus*Stt), a*(Stt - nus*Srr), 2*a*(1+nus)*Srt

Crr, Ctt, Cgrt = lead_coeffs(th, T, K)
print(f"   leading: C_rr = {Crr:.6f}, C_tt = {Ctt:.6f}, C_grt = {Cgrt:.6f}  mm^1/2")

# ---- 3. C_y(90 deg) two independent ways -----------------------------------
Crr90 = lead_coeffs(np.pi/2, T, K)[0]
Cy_disp = (2*K/T)/np.sqrt(2*np.pi) * np.sin(np.pi/4) * (2 - (1+nus)*np.cos(np.pi/4)**2)
print(f"\n3) C_y(90) via displacement formula = {Cy_disp:.6f}; via 2*C_rr(90) = {2*Crr90:.6f} mm^1/2")

# ---- 4. screening interval on theta = 0 ------------------------------------
r_lo = (abs(CZZ_MEAS) / 0.05) ** 2                     # max|principal strain| = |eps_zz| = 0.05
r_hi = J * T / (2*np.pi*(10*sy)**2)                    # sigma_e = 10 sigma_y
print(f"\n4) strain screen r >= {r_lo:.6f}; sigma_e screen r <= {r_hi:.6f}; cap r <= 0.2")
print(f"   -> interval [{min(r_lo, r_hi):.6f}, {min(r_hi, 0.2):.6f}] mm"
      f"  (empty if r_lo > min(r_hi,0.2): {r_lo > min(r_hi, 0.2)})")
for r in (r_lo, r_hi):
    sxx, syy, sxy = sig_K(r, 0.0, K)
    exx, eyy, exy, ezz, se = full_strain(sxx, syy, sxy, T)
    print(f"   r={r:.5f}: sigma_e={se:.2f}, max|princ strain|={max(abs(exx),abs(eyy),abs(ezz)):.5f}")

# ---- verification A: full-law C_zz convergence ------------------------------
print("\nV-A) sqrt(r)*eps_zz(r,0) from FULL law (should -> -0.005 like const*sqrt(r)):")
rs = np.array([0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001])
vals = []
for r in rs:
    sxx, syy, sxy = sig_K(r, 0.0, K)
    exx, eyy, exy, ezz, se = full_strain(sxx, syy, sxy, T)
    vals.append(np.sqrt(r)*ezz)
    print(f"   r={r:7.4f}: sqrt(r) eps_zz = {np.sqrt(r)*ezz: .6f}   sigma_e = {se:7.2f}")
errs = np.array(vals) - CZZ_MEAS
slope = np.polyfit(np.sqrt(rs), errs, 1)[0]
print(f"   error ~ {slope:.5f}*sqrt(r); extrapolated error at r=0: {np.polyfit(np.sqrt(rs), errs, 1)[1]:.2e}")

# ---- verification B: J of the leading field on circles (real units) ---------
def u_modeI(r, th, K, E_, nu_):
    f = 2*K/E_ * np.sqrt(r/(2*np.pi))
    ux = f*np.cos(th/2)*(1 - nu_ + (1+nu_)*np.sin(th/2)**2)
    uy = f*np.sin(th/2)*(2 - (1+nu_)*np.cos(th/2)**2)
    return ux, uy

def J_circle(R, K, T, nus, k=20001):
    a = np.linspace(-np.pi+1e-9, np.pi-1e-9, k)
    x, y = R*np.cos(a), R*np.sin(a)
    sxx, syy, sxy = sig_K(R, a, K)
    exx = (sxx - nus*syy)/T; eyy = (syy - nus*sxx)/T; exy = (1+nus)*sxy/T
    W = 0.5*(sxx*exx + syy*eyy + 2*sxy*exy)
    h = 1e-6*R
    ux1, uy1 = u_modeI(np.hypot(x+h, y), np.arctan2(y, x+h), K, T, nus)
    ux2, uy2 = u_modeI(np.hypot(x-h, y), np.arctan2(y, x-h), K, T, nus)
    dux = (ux1 - ux2)/(2*h)
    duy = (uy1 - uy2)/(2*h)
    Nx = np.cos(a)
    tx, ty = sxx*Nx + sxy*np.sin(a), sxy*Nx + syy*np.sin(a)
    return np.trapezoid((W*Nx - (tx*dux + ty*duy))*R, a)

print("\nV-B) J of the leading field on circles (must equal supplied J = 0.800):")
for R in (0.02, 0.05, 0.1, 0.2):
    print(f"   R = {R:5.3f} mm: J = {J_circle(R, K, T, nus):.6f} N/mm")

