import numpy as np, json, cmath
from scipy import integrate

nu, nu_u, alpha = 0.2, 0.4, 0.89
G, L, sigma = 1.0, 1.0, 1.0
pi = np.pi
sqrt = np.sqrt

print("=" * 70)
print("(a) Loading 1: sigma_n = -H(t), p = 0 on crack faces (permeable crack)")
print("=" * 70)

# Undrained (t=0+) and drained (t->inf) limits are elastic crack problems with
# Poisson's ratio nu_u and nu respectively.  For a Griffith crack under uniform
# face pressure, K_I = p*sqrt(pi*L) independent of elastic constants, and
# V = pi*(1-nu)*p*L^2/G (plane strain).  Hence K is the same at 0+ and inf.
K1_0 = sigma
K1_inf = sigma
V1_0 = pi * (1 - nu_u) * sigma * L**2 / G
V1_inf = pi * (1 - nu) * sigma * L**2 / G

# Short-time: 1D consolidation skin on the crack faces. Undrained pore pressure
# at the face: p* = -(nu_u-nu)*sigma_kk/(alpha*(1-2nu)) with sigma_kk(faces)=-2*sigma.
p_star = 2 * (nu_u - nu) / (alpha * (1 - 2 * nu)) * sigma      # = p_faces_0
A_coeff = 4 * alpha * (1 - 2 * nu) * p_star / ((1 - nu) * sqrt(pi))  # = 8(nu_u-nu)/((1-nu)sqrt(pi))
p_faces_0 = p_star

print(f"K1_0    = {K1_0:.10f}   (x sqrt(pi L), per unit |sigma_n|)")
print(f"K1_inf  = {K1_inf:.10f}")
print(f"V1_0    = {V1_0:.10f}   (x L^2/G)")
print(f"V1_inf  = {V1_inf:.10f}")
print(f"V1_short_coeff A = {A_coeff:.10f}")
print(f"p_faces_0 = {p_faces_0:.10f}")

print()
print("=" * 70)
print("(b) Loading 2: sigma_n = 0, p = H(t) on crack faces")
print("=" * 70)
# t=0+: no mechanical load applied, diffusion layer of zero measure -> K=0, V=0.
# t->inf: steady state p_ss harmonic (p=1 on slit, 0 at inf). The state
# sigma_ij = 0 (identically), e_ij = alpha(1-2nu) p_ss/(2G) delta_ij satisfies
# ALL equations + BCs (strain is compatible since p_ss is harmonic) -> K2_inf=0,
# V2_inf from swelling displacement, tip-closure condition u_y(+-L)=0.
K2_0 = 0.0
K2_inf = 0.0
V2_inf = 2 * alpha * (1 - 2 * nu) * L**2 / (pi * G)
print(f"K2_0    = {K2_0:.10f}")
print(f"K2_inf  = {K2_inf:.10f}")
print(f"V2_inf  = {V2_inf:.10f}")

# ----------------------------------------------------------------------
print()
print("=" * 70)
print("CHECKS")
print("=" * 70)

print("\n[1] Elastic crack via Westergaard Z(z)=tau*(1-z/sqrt(z^2-L^2)),"
      " tau=-1 (pressure loading).")
print("    compare: sigma_yy(x,0) for x>L -> K_I = 1*sqrt(pi L); sigma_xx(face)=-1")
def Zw(z, tau=-1.0, L=1.0):
    return tau * (1.0 - z / np.sqrt(z * z - L * L + 0j))
for eps in [1e-3, 1e-4, 1e-5, 1e-6]:
    z = (L + eps) + 0j
    syy = Zw(z).real
    K_num = syy * sqrt(2 * pi * eps)
    # constant term -tau in Z gives O(sqrt(eps)) bias: extrapolated value:
    K_ext = (K_num - (-1.0) * sqrt(2 * pi * eps)) / sqrt(pi * L)
    print(f"    x-L={eps:.0e}: K_num/sqrt(pi L) = {K_num/sqrt(pi*L):.6f},"
          f" bias-corrected = {K_ext:.8f}  (-> 1.0)")
sx_face = [Zw(x + 1e-9j).real for x in [0.05, 0.3, 0.7, 0.95]]
sy_face = [Zw(x + 1e-9j).real for x in [0.05, 0.3, 0.7, 0.95]]
print(f"    sigma_xx on faces = {['%.6f' % v for v in sx_face]}  (expect -1)")
print(f"    => sigma_kk(faces) = -2  => p_u(face)=2(nu_u-nu)/(alpha(1-2nu)) = "
      f"{2*(nu_u-nu)/(alpha*(1-2*nu)):.8f}")
print("    FAIL MODE: if the Westergaard function or K extraction were wrong,"
      " K_num would not converge to 1 or sigma_xx(face) != -1.")

print("\n[2] Crack volume from COD delta(x)=2(1-nu)*p*sqrt(L^2-x^2)/G:")
for n_, lbl in [(nu_u, "V1_0 (nu_u)"), (nu, "V1_inf (nu)")]:
    cod = lambda x: 2 * (1 - n_) * sigma * sqrt(L**2 - x**2) / G
    V_num, _ = integrate.quad(cod, -L, L)
    print(f"    {lbl}: numeric = {V_num:.10f}, formula pi(1-nu) = {pi*(1-n_):.10f}")
print("    FAIL MODE: mismatch of numeric integral vs closed form would mean wrong COD.")

print("\n[3] Steady loading-2 pore pressure p_ss=(2/pi)Re[arcsin(L/z)]:")
def pss(x, y):
    return (2 / pi) * np.real(cmath.asin(complex(L, 0) / complex(x, y)))
h = 1e-4
lap_pts = []
for (x0, y0) in [(0.5, 0.3), (-0.4, 0.7), (1.5, 0.5), (0.0, 1.2)]:
    lap = (pss(x0 + h, y0) + pss(x0 - h, y0) + pss(x0, y0 + h) + pss(x0, y0 - h)
           - 4 * pss(x0, y0)) / h**2
    lap_pts.append(lap)
print(f"    Laplacian(p_ss) at 4 points = {['%.2e' % v for v in lap_pts]} (expect ~0)")
print(f"    p_ss on slit (x,0+), x=0.3,0.6,0.9: "
      f"{[round(pss(x,1e-9),8) for x in [0.3,0.6,0.9]]} (expect 1)")
print(f"    p_ss at (x,y)=(5,5): {pss(5,5):.2e}; (20,10): {pss(20,10):.2e} (expect ->0)")
print("    FAIL MODE: nonzero Laplacian or wrong boundary values would kill the"
      " sigma=0 steady solution of loading 2.")

print("\n[4] V2_inf: V = 4*c*L^2/pi with c=alpha(1-2nu)/(2G); check via "
      "I=int_0^L x*arccosh(L/x)dx = L^2/2 and u_y(x,0+)=c*(2/pi)*int_x^L arccosh(L/x')dx':")
I_num, _ = integrate.quad(lambda x: x * np.arccosh(np.clip(L / x, 1, None)), 0.0, L)
print(f"    I_num = {I_num:.10f}  vs L^2/2 = 0.5")
def uy(x):
    if x >= L: return 0.0
    inner, _ = integrate.quad(lambda xp: np.arccosh(np.clip(L / xp, 1, None)),
                              x, L, points=[x, L], limit=200)
    return (alpha * (1 - 2 * nu) / (2 * G)) * (2 / pi) * inner
xs = np.linspace(0, L, 401)
uys = np.array([uy(x) for x in xs])
V_num = 2 * 2 * integrate.simpson(uys, x=xs)   # 2 * int_{-L}^{L} u_y dx
print(f"    V2_inf numeric (2*int u_y dx) = {V_num:.8f} vs formula {V2_inf:.8f}")
print("    FAIL MODE: disagreement would mean wrong swelling-displacement construction.")

print("\n[5] Loading-2 short-time skin sign/coeff (not requested, consistency): "
      "dV2 ~ -4*alpha(1-2nu)/((1-nu)sqrt(pi)) * sqrt(ct/L^2):")
print(f"    that coefficient = {-4*alpha*(1-2*nu)/((1-nu)*sqrt(pi)):.8f} (V2 dips <0 then -> +{V2_inf:.4f})")

print("\n[6] skin diffusion p=p**erf(y/(2 sqrt(ct))): check int_0^inf erfc(y/(2 sqrt(ct)))dy = 2 sqrt(ct)/sqrt(pi):")
from scipy.special import erfc
ct = 1e-4
num, _ = integrate.quad(lambda y: erfc(y / (2 * sqrt(ct))), 0, 200 * sqrt(ct))
print(f"    numeric = {num:.8f}, 2*sqrt(ct)/sqrt(pi) = {2*sqrt(ct)/sqrt(pi):.8f}")
print(f"    => dV1 = 4L*alpha(1-2nu)*p0*sqrt(ct)/(G(1-nu)sqrt(pi)) -> A = {A_coeff:.8f}")
print("    FAIL MODE: if the skin constant were wrong, A would not follow from this.")

print("\n[7] Biot consistency: B = 3(nu_u-nu)/(alpha(1-2nu)) = "
      f"{3*(nu_u-nu)/(alpha*(1-2*nu)):.8f}, p_u(face)=2B/3*sigma = {2*(nu_u-nu)/(alpha*(1-2*nu)):.8f}")
print("    Skempton: p_u = -B*sigma_kk/3, sigma_kk(faces) = -2*sigma.")

# ----------------------------------------------------------------------
results = {
    "K1_0": K1_0,
    "K1_inf": K1_inf,
    "V1_0": V1_0,
    "V1_inf": V1_inf,
    "V1_short_coeff": A_coeff,
    "p_faces_0": p_faces_0,
    "K2_0": K2_0,
    "K2_inf": K2_inf,
    "V2_inf": V2_inf,
    "formulas": {
        "K1_0": "sigma*sqrt(pi*L) (undrained elastic limit, nu_u; K of a Griffith crack is independent of elastic constants) -> 1 in units sqrt(pi*L) per unit |sigma_n|",
        "K1_inf": "sigma*sqrt(pi*L) (drained elastic limit, nu) -> 1 in units sqrt(pi*L)",
        "V1_0": "pi*(1-nu_u)*sigma*L^2/G",
        "V1_inf": "pi*(1-nu)*sigma*L^2/G",
        "V1_short_coeff": "A = 8*(nu_u-nu)/((1-nu)*sqrt(pi)), i.e. V ~ V1_0 + A*(L^2/G)*sqrt(c*t/L^2)",
        "p_faces_0": "2*(nu_u-nu)/(alpha*(1-2*nu)) * sigma  (= 2B/3, B=3(nu_u-nu)/(alpha(1-2nu)), since sigma_kk = -2*sigma on the faces)",
        "K2_0": "0 (no mechanical load at t=0+; diffusive boundary layer of vanishing thickness)",
        "K2_inf": "0 (steady state: p_ss=(2/pi)Re[arcsin(L/z)] harmonic and sigma_ij = 0 identically is the exact solution)",
        "V2_inf": "2*alpha*(1-2*nu)*L^2/(pi*G) (swelling displacement u_y=c*(2/pi)*int_x^L arccosh(L/x')dx', c=alpha(1-2nu)/(2G))"
    },
    "reproduce_command": "python3 solve.py"
}
with open("results.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nresults.json written.")
