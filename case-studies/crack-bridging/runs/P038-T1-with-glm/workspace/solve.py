"""Finite slit (half-length a) in infinite plate, plane strain, advancing with
constant bridging traction t on a<|x|<c (=a+L) and intrinsic toughness
K0 = sqrt(E'*Gamma0). Find Sigma and L when delta at x=a (trailing edge of
bridging zone) first reaches delta_c (first connector failure).

Three independent routes:
  A. Tada closed forms (Bao&Suo eqs 45,46) generalized to K_tip = K0 != 0:
       K_tip = sqrt(pi c) [ Sigma - (2t/pi) arccos(a/c) ]            (46)
       delta_t = (4/E') sqrt(c^2-a^2)[Sigma - (2t/pi) arccos(a/c)] + (8 t a/pi E') ln(1+L/a)
              = (4 K0/(E' sqrt(pi c))) sqrt(c^2-a^2) + (8 t a/(pi E')) ln(1+L/a)
  B. Logarithmic-kernel COD quadrature on the finite crack (independent of A):
       delta(x) = (2/(pi E')) Int_{-c}^{c} q(s) G(x,s) ds,
       G = ln|(sqrt(c^2-x^2)+sqrt(c^2-s^2))/(sqrt(c^2-x^2)-sqrt(c^2-s^2))|,
       q = Sigma (whole slit) - t on zones; Sigma from A's K condition.
  C. Spectral dislocation solver from ./knowledge (solve_bridged), linear
     superposition of unit-remote and unit-bridging states to impose K_tip=K0.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
import sys
sys.path.insert(0, 'knowledge')
from crack_bridging import solve_bridged, trailing_edge_theta

# ---------------- constants ----------------
E    = 130e9        # Pa
nu   = 0.28
Ep   = E/(1-nu**2)  # plane strain
a    = 2.40e-3      # m
t    = 160e6        # Pa
dc   = 9.00e-6      # m
G0   = 22.0         # J/m^2
K0   = np.sqrt(Ep*G0)

def Sigma_of_c(c):           # K_tip = K0 condition
    return K0/np.sqrt(np.pi*c) + (2*t/np.pi)*np.arccos(a/c)

def delta_t_A(L):            # route A closed form
    c = a + L
    bracket = K0/np.sqrt(np.pi*c)                      # = Sigma - (2t/pi)arccos
    return (4*bracket/Ep)*np.sqrt(c**2 - a**2) + (8*t*a/(np.pi*Ep))*np.log(1 + L/a)

def delta_t_B(L):            # route B: ln-kernel quadrature, one adaptive call
    c = a + L
    S = Sigma_of_c(c)
    Rx = np.sqrt(c**2 - a**2)
    def G(s):
        Rs = np.sqrt(c**2 - s**2)
        return np.log((Rx + Rs)/abs(Rx - Rs))
    # log singularity at s=a (integrable), G->0 at s=c; split to help quad
    tot = quad(G, a, c, points=[a, 0.5*(a+c)], limit=100,
               epsabs=1e-13, epsrel=1e-12)[0]
    # remote part over whole slit integrates to 4*Sigma*sqrt(c^2-a^2)/Ep (exact)
    return (4*S/Ep)*Rx - (4*t/(np.pi*Ep))*tot

def delta_t_C(L, N=400):     # route C: spectral solver + superposition
    c = a + L
    sr = solve_bridged(L, s=0.0, sigma0=0.0, remote=1.0, N=N, E=Ep, c=c)
    sb = solve_bridged(L, s=0.0, sigma0=1.0, remote=0.0, N=N, E=Ep, c=c)
    rem = (K0 - t*sb['K_tip'])/sr['K_tip']
    tha = trailing_edge_theta(L, c=c)
    dt = rem*sr['delta_at'](np.array([tha]))[0] + t*sb['delta_at'](np.array([tha]))[0]
    return dt, rem, sr['K_tip'], sb['K_tip']

LA = brentq(lambda L: delta_t_A(L) - dc, 1e-6, 0.2, xtol=1e-14, rtol=1e-14)
LB = brentq(lambda L: delta_t_B(L) - dc, 1e-6, 0.2, xtol=1e-14, rtol=1e-14)
LC = brentq(lambda L: delta_t_C(L)[0] - dc, 1e-6, 0.2, xtol=1e-14, rtol=1e-14)

for name, L in (('A', LA), ('B', LB), ('C', LC)):
    c = a + L
    S = Sigma_of_c(c)
    print(f"route {name}: L = {L*1e3:.6f} mm   Sigma = {S/1e6:.5f} MPa   c = {c*1e3:.4f} mm")

# route C internal check of Sigma consistency
dtC, remC, Kr, Kb = delta_t_C(LA)
print(f"  spectral Sigma at L_A: {remC/1e6:.5f} MPa vs analytic {Sigma_of_c(a+LA)/1e6:.5f} MPa"
      f"  (rel diff {abs(remC-Sigma_of_c(a+LA))/Sigma_of_c(a+LA):.2e})")

# ---------------- validations ----------------
print("\n--- validation 1: K0=0 strength state vs knowledge solver free_remote ---")
for L in (1e-3, 3e-3, 6e-3):
    c = a + L
    sol = solve_bridged(L, s=0.0, sigma0=1.0, free_remote=True, N=400, E=Ep, c=c)
    S_ana = (2*t/np.pi)*np.arccos(a/c)
    print(f"  L={L*1e3:.1f} mm: Sigma analytic {S_ana/1e6:.4f} vs spectral {t*sol['remote']/1e6:.4f} MPa")

print("\n--- validation 2: SSB/Dugdale limits of route-A formulas (L << a) ---")
L = 1e-5
print(f"  -K_tip_bridge/sqrt(8t^2 L/pi) = "
      f"{(np.sqrt(np.pi*(a+L))*(2*t/np.pi)*np.arccos(a/(a+L)))/(t*np.sqrt(8*L/np.pi)):.6f} (->1)")
Kb_num = np.sqrt(np.pi*(a+L))*(2*t/np.pi)*np.arccos(a/(a+L))
print(f"  delta_t*Ep/(8 t L) = {delta_t_A(L)*Ep/(8*t*L):.6f} (->1)")

print("\n--- validation 3: eq (47) LSB closed form (K0=0, delta_t=delta_c) ---")
ell = Ep*dc/t
alpha = a/ell
S47 = t*(2/np.pi)*np.arccos(np.exp(-np.pi/(8*alpha)))
L47 = a*(np.exp(np.pi*Ep*dc/(8*t*a)) - 1)
def d47(L):  # route A with K0=0
    c = a+L
    return (8*t*a/(np.pi*Ep))*np.log(1+L/a)
print(f"  alpha = a/ell = {alpha:.4f}; eq(47) Sigma = {S47/1e6:.4f} MPa at L = {L47*1e3:.4f} mm")

print("\n--- validation 4: route B vs route A along the branch ---")
for L in (0.5e-3, 2e-3, LA, 8e-3):
    print(f"  L={L*1e3:8.4f} mm  A: {delta_t_A(L)*1e6:9.5f} um   B: {delta_t_B(L)*1e6:9.5f} um"
          f"   C: {delta_t_C(L)[0]*1e6:9.5f} um")

print("\n--- validation 5: branch monotonicity ---")
Ls = np.linspace(1e-4, 8e-3, 12)
Ss = [Sigma_of_c(a+L)/1e6 for L in Ls]
print("  Sigma(MPa):", " ".join(f"{s:.2f}" for s in Ss))

print("\n--- final ---")
c = a + LA
print(f"L  = {LA*1e3:.4f} mm per end")
print(f"Sigma = {Sigma_of_c(c)/1e6:.3f} MPa")
print(f"check delta_t(A) = {delta_t_A(LA)*1e6:.6f} um; B: {delta_t_B(LA)*1e6:.6f} um; C: {delta_t_C(LA)[0]*1e6:.6f} um")
