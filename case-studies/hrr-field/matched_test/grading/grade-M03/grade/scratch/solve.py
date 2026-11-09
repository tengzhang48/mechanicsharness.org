import numpy as np
from scipy.optimize import brentq

E = 70000.0
nu = 0.30
sy = 10.0
Czz = -0.005
J = 0.8

# Universal mode I angular stress functions (plane, isotropic linear elastic; nu-independent stresses)
def f_stress(theta):
    c = np.cos(theta/2.0); s = np.sin(theta/2.0)
    c3 = np.cos(3*theta/2.0); s3 = np.sin(3*theta/2.0)
    frr = 5.0/4.0*c - 1.0/4.0*c3
    ftt = 3.0/4.0*c + 1.0/4.0*c3
    frt = 1.0/4.0*s + 1.0/4.0*s3
    # Cartesian, for cross-check
    fxx = c*(1 - s*s3)
    fyy = c*(1 + s*s3)
    fxy = s*c*c3
    return frr, ftt, frt, fxx, fyy, fxy

# Effective isotropic moduli of the asymptotic (r->0, <1-sy/se>->1) law
def eff_moduli(T):
    Delta = 1.0/T - 1.0/E
    C_dev = (1.0+nu)/E + 1.5*Delta      # = 1/(2 G_eff)
    C_bulk = (1.0-2*nu)/E               # = 1/(3 K)
    G_eff = 1.0/(2*C_dev)
    K_bulk = 1.0/(3*C_bulk)
    E_eff = 1.0/((2.0/3.0)*C_dev + (1.0/3.0)*C_bulk)
    nu_eff = (3*K_bulk - 2*G_eff)/(2*(3*K_bulk + G_eff))
    return Delta, G_eff, E_eff, nu_eff, C_dev, C_bulk

# A_tr: coefficient in  Czz = -(2K/sqrt(2pi)) * A_tr
def A_tr(T):
    Delta = 1.0/T - 1.0/E
    return nu/E + 0.5*Delta

target = Czz**2 * np.pi / (2.0*J)
print("target A_tr^2 E_eff =", target)

def F(T):
    _, G_eff, E_eff, nu_eff, _, _ = eff_moduli(T)
    return A_tr(T)**2 * E_eff - target

T = brentq(F, 1.0, E-1.0, xtol=1e-12, rtol=1e-15)
Delta, G_eff, E_eff, nu_eff, C_dev, C_bulk = eff_moduli(T)
print("T =", T)
print("A_tr =", A_tr(T), "E_eff =", E_eff, "nu_eff =", nu_eff, "G_eff =", G_eff)

# K from Czz
K = -Czz*np.sqrt(2*np.pi)/(2.0*A_tr(T))
print("K =", K)
print("check J = K^2/E_eff =", K**2/E_eff)

# ---- strains at theta=60 deg via effective law, direct from original compliance form
th = np.deg2rad(60.0)
frr, ftt, frt, fxx, fyy, fxy = f_stress(th)
D = 1.0/T - 1.0/E
# original law, plane stress sigma_zz=0, sigma_kk=frr+ftt (times K/sqrt(2pi r))
s_kk = (frr+ftt)
srr_dev = (2*frr - ftt)/3.0
stt_dev = (2*ftt - frr)/3.0
def strain_from_law(T):
    D = 1.0/T - 1.0/E
    err = (1+nu)/E*frr - nu/E*s_kk + 1.5*D*srr_dev
    ett = (1+nu)/E*ftt - nu/E*s_kk + 1.5*D*stt_dev
    ert_tensor = (1+nu)/E*frt + 1.5*D*frt   # deviatoric shear = frt
    return err, ett, ert_tensor
err, ett, ert = strain_from_law(T)
pref = K/np.sqrt(2*np.pi)
C_rr = pref*err
C_tt = pref*ett
C_gr = pref*2*ert
print("f_rr,f_tt,f_rt =", frr, ftt, frt)
print("C_rr =", C_rr)
print("C_tt =", C_tt)
print("C_gamma_rtheta =", C_gr)

# cross-check with effective E,nu
err2 = (frr - nu_eff*ftt)/E_eff
ett2 = (ftt - nu_eff*frr)/E_eff
ert2 = (1+nu_eff)/E_eff*frt
print("cross strain coeff:", pref*err2, pref*ett2, pref*2*ert2)

# ---- Part 3: C_y(90)
kappa = (3-nu_eff)/(1+nu_eff)
th90 = np.pi/2
C_y = K*kappa/(4*G_eff*np.sqrt(np.pi))
print("kappa_eff =", kappa, "C_y(90) =", C_y)

# independent check: integrate displacement from strain angular functions on theta grid
thg = np.linspace(0, np.pi, 200001)
frrg, fttg, frtg, *_ = f_stress(thg)
# angular strain functions (times r^q with q=-1/2): eps = Ctheta * r^{-1/2}
def ang_strains(T):
    D = 1.0/T - 1.0/E
    sk = frrg+fttg
    er = (1+nu)/E*frrg - nu/E*sk + 1.5*D*(2*frrg-fttg)/3.0
    et = (1+nu)/E*fttg - nu/E*sk + 1.5*D*(2*fttg-frrg)/3.0
    ert = ((1+nu)/E + 1.5*D)*frtg
    return er, et, ert
cr, ct, crt = ang_strains(T)
# u_r = r^{1/2} U(θ), U_r = 2*C_rr(θ) because eps_rr = d/dr[r^1/2 U] = 1/2 r^{-1/2}U
Ur = 2*cr
# eps_tt = (u_r + du_t/dθ)/r  => U_t' = eps_tt_ang - U_r
integ = np.concatenate([[0.0], np.cumsum(0.5*(np.diff(ct-Ur)*(thg[1]-thg[0]))*0 )])  # placeholder
# do cumulative trapezoid properly
dut = ct - Ur
Ut = np.concatenate([[0.0], np.cumsum(0.5*(dut[1:]+dut[:-1])*np.diff(thg))])
i90 = np.argmin(np.abs(thg-th90))
pref2 = K/np.sqrt(2*np.pi)
# u_y = u_r sinθ + u_θ cosθ; at θ=90 => u_y = u_r*1 = r^1/2 * pref2 * Ur  (since U_r=-Ur? careful)
uy_factor = pref2*Ur[i90]*np.sin(thg[i90]) + pref2*Ut[i90]*np.cos(thg[i90])
print("integrated C_y(90) =", uy_factor)

# ---- Part 4 screens on theta=0
# sigma_e(0) = K/sqrt(2 pi r); require >= 10 sy  => r <= r1
r1 = K**2/(2*np.pi*(10*sy)**2)
# strains at theta=0
frr0, ftt0, frt0, fxx0, fyy0, fxy0 = f_stress(0.0)
D = 1.0/T - 1.0/E
sk0 = 2*frr0  # frr=ftt=1, s_kk=2
# in-plane strain (== principal in-plane), out-of-plane strain
e_in_ang = (1+nu)/E*frr0 - nu/E*sk0 + 1.5*D*(2*frr0-ftt0)/3.0
e_zz_ang = -nu/E*sk0 + 1.5*D*(-sk0/3.0)
pref0 = K/np.sqrt(2*np.pi)
# strain(r,0) = pref0/sqrt(r) * ang
def e_in(r): return pref0/np.sqrt(r)*e_in_ang
def e_zz(r): return pref0/np.sqrt(r)*e_zz_ang
print("e_in_ang, e_zz_ang =", e_in_ang, e_zz_ang)
# require max(|e_in|,|e_zz|) <= 0.05
r_in = (pref0*abs(e_in_ang)/0.05)**2
r_zz = (pref0*abs(e_zz_ang)/0.05)**2
r2 = max(r_in, r_zz)
print("r1 (sigma_e>=10sy) =", r1)
print("r_in for 0.05 =", r_in, "r_zz =", r_zz, "r2 =", r2)
r_hi = min(r1, 0.200)
print("interval = [", r2, ",", r_hi, "]  empty?", r2 > r_hi)
print("principal strain values at r2:", e_in(r2), e_zz(r2))
print("sigma_e at r_hi / sy =", (K/np.sqrt(2*np.pi*r_hi))/sy if r_hi>0 else None)
