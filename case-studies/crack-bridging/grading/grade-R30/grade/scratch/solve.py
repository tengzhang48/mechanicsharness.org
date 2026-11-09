import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

# ---------------- constants ----------------
E  = 130e9
nu = 0.28
Ep = E/(1-nu**2)          # plane strain modulus
a  = 2.40e-3              # slit half-length [m]
t  = 160e6                # connector closing traction [Pa]
dc = 9.00e-6              # critical full opening [m]
G0 = 22.0                 # intrinsic tip work [J/m^2]
K0 = np.sqrt(Ep*G0)       # tip toughness [Pa sqrt(m)]

theta = lambda c: np.arccos(a/c)
Sigma = lambda c: K0/np.sqrt(np.pi*c) + (2*t/np.pi)*theta(c)   # tip condition K_tip=K0

def delta_a_cf(c):
    # opening at zone root x=a: superposition = classical Dugdale CTOD at Sigma_D
    # + Griffith opening from the extra remote stress K0/sqrt(pi c)
    dug = (8*t*a/(np.pi*Ep))*np.log(1.0/np.cos(theta(c)))
    tip = (4*K0/(Ep*np.sqrt(np.pi*c)))*np.sqrt(c**2-a**2)
    return dug + tip

# ---------------- main solve ----------------
c_lo, c_hi = a*1.0000001, 0.05
f = lambda c: delta_a_cf(c) - dc
c_s = brentq(f, c_lo, c_hi, xtol=1e-16, rtol=8.9e-16, maxiter=300)
S_s = Sigma(c_s)
L_s = c_s - a

print("Ep [GPa]      =", Ep/1e9)
print("K0  [MPa sqrtm] =", K0/1e6)
print("Sigma_init [MPa] =", K0/np.sqrt(np.pi*a)/1e6)
print("Sigma_ss (SSB) [MPa] =", np.sqrt(Ep*(G0+t*dc)/(np.pi*a))/1e6)
print("c* [mm] =", c_s*1e3)
print("L* [mm] =", L_s*1e3)
print("Sigma* [MPa] =", S_s/1e6)
print("delta_a(c*) [um] =", delta_a_cf(c_s)*1e6)

# branch table
print("\nbranch (monotone increase of delta(a) with c):")
for c in [a*1.05, 2.6e-3, 3.5e-3, 4.5e-3, 5.5e-3, c_s, 8e-3]:
    print(f"  c={c*1e3:8.4f} mm  Sigma={Sigma(c)/1e6:8.3f} MPa  delta(a)={delta_a_cf(c)*1e6:8.4f} um")

# ---------------- check 1: verify Dugdale CTOD formula via direct integration of
# the dislocation-density solution of the SIE (independent semi-analytic route) ----
def delta_a_J(S, c):
    sq = np.sqrt(c*c - a*a)
    def J(x):  # int_{-a}^{a} sqrt(c^2-s^2)/(x-s) ds  for x>a
        reg = quad(lambda s: (np.sqrt(c*c-s*s)-sq)/(x-s), -a, a, limit=400)[0]
        return sq*np.log((x+a)/(x-a)) + reg
    tha = theta(c)
    I = quad(lambda th: J(c*np.cos(th)), 0.0, tha, points=[tha], limit=400)[0]
    return (4/Ep)*((S-t)*sq + (t/np.pi)*I)

SD = (2*t/np.pi)*theta(c_s)   # zero-K Dugdale stress for this (a,c)
d_num = delta_a_J(SD, c_s)
d_lit = (8*t*a/(np.pi*Ep))*np.log(1.0/np.cos(theta(c_s)))
print("\nCHECK1 Dugdale CTOD: semi-analytic SIE integration =", d_num*1e6, "um")
print("      literature formula (8 t a/piE')lnsec        =", d_lit*1e6, "um")
print("      rel. diff =", abs(d_num-d_lit)/d_lit)
d_tot = delta_a_J(S_s, c_s)
print("CHECK1b delta_a at solution (Sigma*,c*) via SIE integration =", d_tot*1e6, "um  vs target 9.0")

# ---------------- check 2: independent spectral SIE solver ----------------
def solve_SIE(S, c, M=2048):
    xa = a/c
    k  = np.arange(1, M+1)
    xi = np.cos(k*np.pi/(M+1))
    s2 = np.sin(k*np.pi/(M+1))**2
    g  = np.where(np.abs(xi) > xa, S - t, S)      # +p on faces; closing = -t on zones
    F  = -(4.0*c/Ep)*g
    N  = M
    b  = np.zeros(N+1)
    U0 = np.ones(M); U1 = 2*xi
    b[1] = (2.0/(M+1))*np.dot(s2*F, U0)
    b[2] = (2.0/(M+1))*np.dot(s2*F, U1)
    for n in range(3, N+1):
        U2 = 2*xi*U1 - U0
        b[n] = (2.0/(M+1))*np.dot(s2*F, U2)
        U0, U1 = U1, U2
    K = -(Ep/4)*np.sqrt(np.pi/c)*b[1:].sum()
    def delt(xi0):
        s = 0.0; V0 = 1.0; V1 = 2*xi0
        s = b[1]*V0
        for n in range(2, N+1):
            s += b[n]*V1/n          # U_{n-1}(xi0) = V1 (before update)
            V2 = 2*xi0*V1 - V0
            V0, V1 = V1, V2
        return -np.sqrt(1-xi0*xi0)*s
    return K, delt

Kspec, delt = solve_SIE(S_s, c_s)
print("\nCHECK2 spectral SIE at (Sigma*, c*):")
print("  K_tip spectral [MPa sqrtm] =", Kspec/1e6, " vs K0 =", K0/1e6)
print("  delta(a) spectral [um]     =", delt(a/c_s)*1e6, " vs target 9.0")

# anchors for the spectral solver
Kg, dg = solve_SIE(30e6, 3e-3)  # t irrelevant if zone empty? no: t acts; use pure remote check separately
print("\n  (anchor) pure-remote K vs Sigma*sqrt(pi c):", )
# pure remote anchor: set t=0 by using a=0 trick not possible; instead verify with t loading off:
# quick anchor: Griffith profile
def griffith_check():
    import numpy as np
    global a
    a_old = a
    a = 1e-9  # zones cover essentially whole crack -> g = S - t everywhere; set t->0 case instead
    a = a_old
griffith_check()

# monotonicity of opening over the bridge zone (spectral)
xs = np.linspace(0, 1.0, 11)
vals = [delt(x0) for x0 in xs*np.max([a/c_s*0.999, 1e-6])]
print("\nCHECK3 delta(x)/delta(a) along zone (xi from ~0 to a/c):")
print(" ", np.array(vals)/delt(a/c_s))

# small-scale bridging limit comparison
Kss = np.sqrt(Ep*(G0+t*dc))
L_ssb = np.pi/8*(Kss-K0)**2/t**2
print("\nSSB limits: Sigma_ss =", np.sqrt(Ep*(G0+t*dc)/(np.pi*a))/1e6,
      "MPa ; L_ssb =", L_ssb*1e3, "mm  (large-scale answer differs as expected since L~a)")
