"""Infer T and near-tip coefficients for the bilinear total-strain (deformation-theory) law.

Reduction: as r->0, sigma_e ~ r^-1/2 -> inf, so <1 - sigma_y/sigma_e> -> 1 and
    eps_ij = (1+nu*)/E* sigma_ij - (nu*/E*) sigma_kk delta_ij  (plane-stress 2D law, eps_zz free)
with E* = T exactly (since eps_kk = (1-2nu)/E sigma_kk is purely elastic) and
    nu* = 1/2 - (T/(2E))(1-2nu).
The leading field is the plane-stress Williams mode-I field with (E*, nu*) = (T, nu*), and
J = K^2/T (plane stress).  C_zz(0) = -(2 nu*/T) K/sqrt(2 pi).
"""
import numpy as np
from scipy.optimize import brentq

E, nu, sy = 70000.0, 0.30, 10.0
J = 0.800                    # N/mm
Czz0 = -0.005000             # mm^1/2, measured

def nu_star(T):  return 0.5 - (T / (2 * E)) * (1 - 2 * nu)
def K_of(T):     return np.sqrt(J * T)

def eq(T):  # leading-field C_zz(0) minus measured
    K = K_of(T)
    return -(2 * nu_star(T) / T) * K / np.sqrt(2 * np.pi) - Czz0

T = brentq(eq, 1.0, E - 1.0, xtol=1e-12, rtol=8.9e-16)
K = K_of(T)
nus = nu_star(T)
print(f"T   = {T:.6f} MPa   (E* = T, nu* = {nus:.8f})")
print(f"K   = {K:.6f} MPa*mm^1/2   from J = K^2/T = {K**2/T:.6f}")
print(f"check: sigma_y/E + (sigma-sigma_y)/T bilinear;  T < E: {T < E}")

# -------------------------------- Williams plane-stress leading field, angular functions
def f_ang(th):
    """Angular stress functions: sigma_ij = K/sqrt(2 pi r) * f_ij(th), polar components."""
    c2, s2 = np.cos(th / 2), np.sin(th / 2)
    sxx = c2 * (1 - s2 * np.sin(3 * th / 2))
    syy = c2 * (1 + s2 * np.sin(3 * th / 2))
    sxy = c2 * s2 * np.cos(3 * th / 2)
    c, s = np.cos(th), np.sin(th)
    srr = sxx * c**2 + syy * s**2 + 2 * sxy * s * c
    stt = sxx * s**2 + syy * c**2 - 2 * sxy * s * c
    srt = (syy - sxx) * s * c + sxy * (c**2 - s**2)
    return srr, stt, srt

def coeffs(th):
    """C_rr, C_tt, C_gamma_rt = sqrt(r)*(eps_rr, eps_tt, 2 eps_rt) of the leading field."""
    srr, stt, srt = f_ang(np.asarray(th, float))
    a = K / (T * np.sqrt(2 * np.pi))
    return a * (srr - nus * stt), a * (stt - nus * srr), a * 2 * (1 + nus) * srt

th60 = np.pi / 3
Crr, Ctt, Cgrt = (float(v) for v in coeffs(th60))
print(f"\ntheta = 60 deg: C_rr = {Crr:.6f}, C_tt = {Ctt:.6f}, C_gamma_rt = {Cgrt:.6f} mm^1/2")

# C_y(90°): u_y = u_r at theta = 90°, and u_r = sqrt(r)*U_r with U_r = 2 e_rr_ang (q = -1/2)
Cy = 2 * float(coeffs(np.pi / 2)[0])
print(f"C_y(90°) = {Cy:.6f} mm^1/2   [kinematics: u_y(90°) = u_r(90°), U_r = 2 e_rr^ang]")

# -------------------------------- screens on theta = 0
# sigma_e(r,0) = K/sqrt(2 pi r) (plane stress, sigma_xx = sigma_yy);  principal strains:
# eps_1 = eps_2 = (1-nu*)sigma/T, eps_3 = -2 nu* sigma/T, sigma = K/sqrt(2 pi r)
sig = lambda r: K / np.sqrt(2 * np.pi * r)
e_pr = lambda r: np.array([(1 - nus) * sig(r) / T, (1 - nus) * sig(r) / T, -2 * nus * sig(r) / T])
r_strain = (abs(Czz0) / 0.0500) ** 2          # |eps_zz| = |Czz0|/sqrt(r) <= 0.05
r_yeild = (K / (np.sqrt(2 * np.pi) * 10 * sy)) ** 2   # sigma_e >= 10 sigma_y
r_lo = max(r_strain, 0.0)
r_hi = min(r_yeild, 0.200)
print(f"\nscreens: r >= {r_strain:.6f} (strain), r <= {r_yeild:.6f} (sigma_e), r <= 0.200")
print(f"interval = [{r_lo:.6f}, {r_hi:.6f}] mm, nonempty = {r_lo <= r_hi}")
print(f"  at r_lo: max|princ strain| = {np.max(np.abs(e_pr(r_lo))):.5f}, sigma_e = {sig(r_lo):.2f}")
print(f"  at r_hi: max|princ strain| = {np.max(np.abs(e_pr(r_hi))):.5f}, sigma_e = {sig(r_hi):.2f}")

# -------------------------------- verification 1: full given law vs effective linear law near tip
def eps_total(sig3, T_):
    """Full given law for a 3-component stress tuple (sxx, syy, sxy), plane stress szz=0."""
    sxx, syy, sxy = sig3
    skk = sxx + syy
    s = np.array([sxx - skk / 3, syy - skk / 3, -skk / 3, sxy])  # xx,yy,zz,xy deviator
    se = np.sqrt(1.5 * (s[0]**2 + s[1]**2 + s[2]**2 + 2 * s[3]**2))
    mac = max(1.0 - sy / se, 0.0) if se > 0 else 0.0
    def eps(i): return (1 + nu) / E * sig3[i % 3 if i < 3 else 2]  # unused
    e_xx = (1 + nu) / E * sxx - nu / E * skk + 1.5 * (1 / T_ - 1 / E) * mac * s[0]
    e_yy = (1 + nu) / E * syy - nu / E * skk + 1.5 * (1 / T_ - 1 / E) * mac * s[1]
    e_zz = -nu / E * skk + 1.5 * (1 / T_ - 1 / E) * mac * s[2]
    e_xy = (1 + nu) / E * sxy + 1.5 * (1 / T_ - 1 / E) * mac * s[3]
    return np.array([e_xx, e_yy, e_xy, e_zz])

print("\nfull-law vs effective-linear-law relative gap at theta=60° (should be ~ sigma_y/sigma_e):")
for r in (0.010, 0.030, 0.061):
    th = th60
    c2, s2 = np.cos(th / 2), np.sin(th / 2)
    f = K / np.sqrt(2 * np.pi * r)
    sig3 = (f * c2 * (1 - s2 * np.sin(3 * th / 2)), f * c2 * (1 + s2 * np.sin(3 * th / 2)),
            f * c2 * s2 * np.cos(3 * th / 2))
    skk = sig3[0] + sig3[1]
    se = np.sqrt(sig3[0]**2 + sig3[1]**2 - sig3[0] * sig3[1] + 3 * sig3[2]**2)
    e_tot = eps_total(sig3, T)
    e_lin = np.array([(sig3[0] - nus * sig3[1]) / T, (sig3[1] - nus * sig3[0]) / T,
                      (1 + nus) * sig3[2] / T, -nus / T * skk])
    gap = np.max(np.abs(e_tot - e_lin)) / np.max(np.abs(e_lin))
    print(f"  r = {r:.3f}: sigma_e = {se:8.2f} MPa, sigma_y/sigma_e = {sy/se:.4f}, "
          f"max strain gap = {gap:.4f}")

# uniaxial sanity: bilinear stress-strain
s = 2 * sy
eps_u = s / E + 1.5 * (1 / T - 1 / E) * (1 - sy / s) * (2 * s / 3)
print(f"\nuniaxial at 2*sigma_y: eps = {eps_u:.6f} vs sigma_y/E + sigma_y/T = "
      f"{sy/E + sy/T:.6f} (bilinear check)")
