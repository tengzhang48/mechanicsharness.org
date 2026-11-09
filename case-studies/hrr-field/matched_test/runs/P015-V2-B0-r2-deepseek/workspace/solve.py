import numpy as np
from scipy.optimize import brentq

E = 70000.0
nu = 0.30
sy = 10.0
J = 0.80
Czz = -0.00500

# Near-tip fully-plastic asymptotic law is linear with uniaxial modulus E_eff = T
# (derive/check below). Plane stress: J = Ktip^2 / E_eff = Ktip^2 / T.
# At theta=0: sigma_e = sigma_xx = sigma_yy = K/sqrt(2 pi r).
# eps_zz = -(2 nu/E) sig - (1/T - 1/E)(sig - sy)  -> asymptotic coeff -P sig, P = 2nu/E + 1/T - 1/E.
def P_of(T):
    return 2*nu/E + 1.0/T - 1.0/E

def Czz_of_T(T):
    K = np.sqrt(J*T)
    return -P_of(T)*K/np.sqrt(2*np.pi)

# solve Czz_of_T(T) = Czz
T = brentq(lambda T: Czz_of_T(T) - Czz, 1.0, E-1e-6, xtol=1e-12, rtol=1e-15)
K = np.sqrt(J*T)

D = (1+nu)/E + 1.5*(1/T - 1/E)      # deviatoric compliance (1/2G_eff)
B = (1-2*nu)/E                       # spherical compliance-ish (1-2nu)/E
E_eff = 3.0/(2*D + B)
nu_eff = (D - B)/(2*D + B)
G_eff = E_eff/(2*(1+nu_eff))
kappa = (3-nu_eff)/(1+nu_eff)

print("T solved        =", T)
print("E_eff           =", E_eff, " (should equal T)")
print("nu_eff          =", nu_eff)
print("K_tip           =", K)
print("P               =", P_of(T))
print("Czz check       =", Czz_of_T(T))

# ---- direct numerical check that J from the effective-linear mode-I field equals K^2/T
def stress_cart(x, y, K):
    r = np.hypot(x, y); th = np.arctan2(y, x)
    c, s = np.cos(th/2), np.sin(th/2)
    fxx = c*(1 - s*np.sin(1.5*th))
    fyy = c*(1 + s*np.sin(1.5*th))
    fxy = s*c*np.cos(1.5*th)
    rr = np.sqrt(x*x+y*y)
    return [K/np.sqrt(2*np.pi*rr)*f for f in (fxx, fyy, fxy)]

def strain_cart(x, y, K):
    sxx, syy, sxy = stress_cart(x, y, K)
    exx = (sxx - nu_eff*syy)/E_eff
    eyy = (syy - nu_eff*sxx)/E_eff
    exy = (1+nu_eff)*sxy/E_eff
    return exx, eyy, exy

def u_cart(x, y, K):
    r = np.hypot(x, y); th = np.arctan2(y, x)
    c, s = np.cos(th/2), np.sin(th/2)
    fac = K/(2*G_eff)*np.sqrt(r/(2*np.pi))
    ux = fac*c*(kappa - 1 + 2*s*s)
    uy = fac*s*(kappa + 1 - 2*c*c)
    return ux, uy

def J_circle(R, K, n=40001):
    a = np.linspace(-np.pi, np.pi, n)
    x, y = R*np.cos(a), R*np.sin(a)
    sxx, syy, sxy = stress_cart(x, y, K)
    exx, eyy, exy = strain_cart(x, y, K)
    W = 0.5*(sxx*exx + syy*eyy + 2*sxy*exy)
    h = 1e-5
    uxp, uyp = u_cart(x+h*R, y, K)
    uxm, uym = u_cart(x-h*R, y, K)
    dux = (uxp-uxm)/(2*h*R)
    duy = (uyp-uym)/(2*h*R)
    Nx, Ny = np.cos(a), np.sin(a)
    tx = sxx*Nx + sxy*Ny
    ty = sxy*Nx + syy*Ny
    return np.trapezoid((W*Nx - (tx*dux + ty*duy))*R, a)

for R in (0.1, 1.0, 10.0):
    print(f"J_num(R={R}) = {J_circle(R,K):.8f}   K^2/T = {K*K/T:.8f}")

# ---- Part 2: theta = 60 deg
th = np.radians(60.0)
c, s = np.cos(th/2), np.sin(th/2)
f_rr = c*(1 + s*s)
f_tt = c*(1 - s*s)
f_rt = s*c*c
C_rr = K/(E_eff*np.sqrt(2*np.pi))*(f_rr - nu_eff*f_tt)
C_tt = K/(E_eff*np.sqrt(2*np.pi))*(f_tt - nu_eff*f_rr)
C_grth = 2*(1+nu_eff)/E_eff*(K/np.sqrt(2*np.pi))*f_rt
print("\nPart2 theta=60deg:")
print(" f_rr,f_tt,f_rt =", f_rr, f_tt, f_rt)
print(" C_rr   =", C_rr)
print(" C_tt   =", C_tt)
print(" C_grth =", C_grth)

# ---- Part 3: C_y(90) = lim u_y/sqrt(r)
ux90, uy90 = u_cart(0.0, 1.0, K)   # r=1 at theta=90
print("\nPart3 C_y(90) =", uy90)

# ---- Part 4: interval on theta=0
r_a = K*K/(2*np.pi*(10*sy)**2)      # sigma_e >= 10 sy
# effective (leading) strains
emax_eff = max(1-nu_eff, 2*nu_eff)
Cs_eff = 0.05*T/emax_eff
r_b_eff = K*K/(2*np.pi*Cs_eff**2)
# full-law strains for comparison: given sigma_xx
def princ_strains_full(sig):
    exx = (1-nu)/E*sig + 0.5*(1/T-1/E)*(sig-sy)
    ezz = -(2*nu/E)*sig - (1/T-1/E)*(sig-sy)
    return exx, exx, ezz
# solve sigma_xx s.t. max|principal| = 0.05 (sigma>>sy)
def maxstrain_full(sig):
    exx, eyy, ezz = princ_strains_full(sig)
    return max(abs(exx), abs(eyy), abs(ezz))
sigB = brentq(lambda s: maxstrain_full(s)-0.05, 100.0, 1e9)
r_b_full = K*K/(2*np.pi*sigB**2)
print("\nPart4:")
print(" r_a (sigma_e>=100)          =", r_a)
print(" r_b effective leading       =", r_b_eff, " (Cs_eff=%.6f)" % Cs_eff)
print(" r_b full-law                =", r_b_full, " (sigB=%.4f)" % sigB)
lo = r_b_full
hi = min(r_a, 0.200)
print(" interval [%g, %g]  nonempty=%s" % (lo, hi, lo <= hi))
print(" principal strain at r=0.2 (full law, sigma=%.4f):" % (K/np.sqrt(2*np.pi*0.2)), princ_strains_full(K/np.sqrt(2*np.pi*0.2)))
print(" sigma_e at r=0.2 =", K/np.sqrt(2*np.pi*0.2))
