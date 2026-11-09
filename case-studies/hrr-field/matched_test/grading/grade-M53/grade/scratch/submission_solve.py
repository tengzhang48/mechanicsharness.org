"""Mode III (anti-plane) power-law crack-tip field, amplitude fixed by J = 0.120 N/mm.

Model: pure power law gamma = (gamma0/tau0)(tau/tau0)^(n-1) tau (deformation theory; proportional
loading inside the near-tip zone makes it exact for the leading term). Separable stress function
psi = a r^m g(theta), m = n/(n+1), i.e. stress exponent p = 1/(n+1) fixed by a finite nonzero J.

  tau_rz = a r^(m-1) g'(th),  tau_thz = -a m r^(m-1) g(th),  rho^2 = g'^2 + m^2 g^2

Compatibility (curl-free gamma) gives the angular ODE
  (rho^(n-1) g')' + m alpha rho^(n-1) g = 0,  alpha = 1/(n+1),
with g(0) = -1 (normalization), g'(0) = 0 (mode III antisymmetry: tau_rz(0)=0), g(pi) = 0
(traction-free faces: tau_thz(+-pi) = 0).

Displacement: w = Cw r^alpha F(th), F = rho^(n-1) g', Cw = (gamma0/alpha/tau0^n) a^n.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ------------------------------------------------------------------ material / loads
tau0 = 80.0      # MPa
gamma0 = 0.002
n = 3.5
J_given = 0.120  # N/mm

alpha = 1.0 / (n + 1.0)
m = n / (n + 1.0)
p = alpha  # stress exponent tau ~ r^-p


def g_ode(mval):
    def rhs(th, y):
        g1, g2 = y
        r2 = g2 * g2 + mval * mval * g1 * g1
        num = g1 * (mval * alpha * r2 + (n - 1.0) * mval * mval * g2 * g2)
        den = r2 + (n - 1.0) * g2 * g2
        return [g2, -num / den]
    return rhs


def integrate_g(mval, thmax=np.pi, rtol=1e-12, atol=1e-14):
    sol = solve_ivp(g_ode(mval), (0.0, thmax), [-1.0, 0.0],
                    dense_output=True, rtol=rtol, atol=atol, max_step=0.01)
    return sol


def shape(mval=m):
    """Return callables g(th), gp(th), rho(th), F(th), Fp(th) on [0, pi]."""
    sol = integrate_g(mval)
    def g(th):   return sol.sol(th)[0]
    def gp(th):  return sol.sol(th)[1]
    def rho(th): return np.hypot(sol.sol(th)[1], mval * sol.sol(th)[0])
    def F(th):   return rho(th) ** (n - 1.0) * gp(th)
    def Fp(th):  return -mval * alpha * rho(th) ** (n - 1.0) * g(th)
    return sol, g, gp, rho, F, Fp


# ------------------------------------------------------------------ J for unit amplitude
def C_J(mval=m):
    """J = C_J * a^(n+1).  Integrand in Cartesian-safe form (see notes)."""
    sol, g, gp, rho, F, Fp = shape(mval)
    x, wq = np.polynomial.legendre.leggauss(400)
    th = 0.5 * np.pi * (x + 1.0)                     # [0, pi]
    wth = 0.5 * np.pi * wq
    gv, gpv, rv = g(th), gp(th), rho(th)
    integrand = (n / (n + 1.0) * rv ** (n + 1.0) * np.cos(th)
                 - rv ** (n - 1.0) * gpv ** 2 * np.cos(th)
                 - mval * gv * gpv * rv ** (n - 1.0) * np.sin(th))
    return 2.0 * np.sum(integrand * wth) * gamma0 / tau0 ** n


# ------------------------------------------------------------------ solve
print("=" * 72)
print("1) eigen-structure checks")
# face residual with m fixed by J-consistency
sol = integrate_g(m)
print(f"   face residual g(pi) with m = n/(n+1) = {m:.10f}:  {sol.y[0, -1]:+.3e}  (need ~0)")
# eigenvalue shot: find m such that g(pi) = 0
f_shot = lambda mv: integrate_g(mv, rtol=1e-11, atol=1e-13).y[0, -1]
m_shot = brentq(f_shot, 0.6, 0.9, xtol=1e-12)
print(f"   eigenvalue from shooting g(pi)=0:           m = {m_shot:.12f}  vs n/(n+1) = {m:.12f}")

# elastic limit n=1: g should be -cos(theta/2), and the J integral I must equal pi
def elastic_check():
    global n, alpha, m
    n_s, a_s, m_s = n, alpha, m
    n, alpha, m = 1.0, 0.5, 0.5
    sol1 = integrate_g(0.5)
    th = np.linspace(0, np.pi, 2001)
    err = np.max(np.abs(sol1.sol(th)[0] + np.cos(th / 2)))
    CJ1 = C_J(0.5)
    I1 = CJ1 * tau0 ** n / gamma0  # = integral value
    # with tau_r = a r^-1/2 g'/2 = (K/sqrt(2 pi)) r^-1/2 sin(th/2)  =>  a = K sqrt(2/pi)
    # elastic check: J = C_J a^2  must equal  K^2/(2 mu),  mu = tau0/gamma0
    a1 = np.sqrt(2 / np.pi)                       # K = 1
    lhs = CJ1 * a1 ** 2; rhs = 1.0 ** 2 / (2 * tau0 / gamma0)
    print(f"   n=1 limit: max|g + cos(th/2)| = {err:.3e};  integral I(n=1) = {I1:.12f} = pi/4 = {np.pi/4:.12f}")
    print(f"   n=1 amplitude check: C_J a^2 = {lhs:.12f}  vs  K^2/(2 mu) = {rhs:.12f}   (rel diff {abs(lhs/rhs-1):.2e})")
    n, alpha, m = n_s, a_s, m_s
elastic_check()

print("=" * 72)
print("2) amplitude from J")
CJ = C_J(m)
a = (J_given / CJ) ** (1.0 / (n + 1.0))
print(f"   C_J = {CJ:.8e} (units: with this normalization J = C_J a^(n+1))")
print(f"   a = {a:.8f} MPa*mm^p   [tau = a r^-p rho]")

# ------------------------------------------------------------------ requested quantities
sol, g, gp, rho, F, Fp = shape(m)
Cw = gamma0 / alpha / tau0 ** n * a ** n

def w_disp(r, th):
    return Cw * r ** alpha * F(th)

def tau_r(r, th):   return a * r ** (m - 1.0) * gp(th)
def tau_th(r, th):  return -a * m * r ** (m - 1.0) * g(th)
def gamma_mag(r, th): return gamma0 * (a * r ** (m - 1.0) * rho(th) / tau0) ** n

print("=" * 72)
print("3) answers")
for r in (0.0200, 0.0800):
    dw = 2.0 * w_disp(r, np.pi)
    print(f"   Delta_w(r = {r:.4f} mm) = {dw:.6e} mm")

r2, th2 = 0.0200, np.deg2rad(60.0)
print(f"   tau_rz(r=0.02, 60deg)  = {tau_r(r2, th2):.6f} MPa")
print(f"   tau_thz(r=0.02, 60deg) = {tau_th(r2, th2):.6f} MPa")

tau_star = tau0 * (0.0100 / gamma0) ** (1.0 / n)
print(f"   |gamma|=0.01 level: tau* = {tau_star:.6f} MPa")
for name, th in (("0 deg", 0.0), ("90 deg", np.pi / 2), ("180 deg (upper face)", np.pi)):
    rstar = (a * rho(th) / tau_star) ** (n + 1.0)
    print(f"   r(|gamma|=0.01, theta={name:20s}) = {rstar:.6f} mm")

# ------------------------------------------------------------------ independent verifications
print("=" * 72)
print("4) verifications")

# (a) J by brute-force Cartesian integration at two radii: FD for w_,x from Cartesian points
#     (w is odd in theta; dense ODE solution is only valid on [0, pi], so integrate theta in [0, pi]
#      and double — the integrand W cos(th) - tau_r w_x is even in theta)
def J_numeric(rcont, h=1e-6, npts=200001):
    th = np.linspace(0.0, np.pi, npts)
    def w_cart(x, y):
        rr = np.hypot(x, y)
        tt = np.arctan2(y, x)
        return w_disp(rr, tt)
    x, y = rcont * np.cos(th), rcont * np.sin(th)
    wxp = w_cart(x + h, y); wxm = w_cart(x - h, y)
    wx = (wxp - wxm) / (2 * h)
    tz = tau_r(rcont, th)                       # traction on circle: sigma_xz Nx + sigma_yz Ny = tau_rz
    Wv = n / (n + 1.0) * gamma0 / tau0 ** n * (a * rcont ** (m - 1.0) * rho(th)) ** (n + 1.0)
    integrand = (Wv * np.cos(th) - tz * wx) * rcont
    return 2.0 * np.trapezoid(integrand, th)

for rc in (0.0200, 0.0800, 0.150):
    print(f"   J on circle r={rc:.3f} (Cartesian FD): {J_numeric(rc):.6f} N/mm   (given 0.120)")

# (b) equilibrium + compatibility residuals by FD in polar (theta interior of [0, pi])
def residuals(r=0.03, nth=4001):
    th = np.linspace(1e-4, np.pi - 1e-4, nth)
    hth = th[1] - th[0]; hr = r * 1e-6
    tr = lambda rr, tt: a * rr ** (m - 1.0) * gp(tt)
    tt_ = lambda rr, tt: -a * m * rr ** (m - 1.0) * g(tt)
    eq = (tr(r + hr, th) * (r + hr) - tr(r - hr, th) * (r - hr)) / (2 * hr * r) \
         + (tt_(r, th + hth) - tt_(r, th - hth)) / (2 * hth * r)
    # compatibility: (1/r) d(r g_th)/dr - d(g_r)/dth = 0 ;  strains ~ r^((m-1)n) = r^(alpha-1)
    gr = lambda rr, tt: gamma0 / tau0 ** n * a ** n * rr ** ((m - 1.0) * n) * rho(tt) ** (n - 1.0) * gp(tt)
    gth = lambda rr, tt: -gamma0 / tau0 ** n * a ** n * rr ** ((m - 1.0) * n) * rho(tt) ** (n - 1.0) * m * g(tt)
    comp = (gth(r + hr, th) * (r + hr) - gth(r - hr, th) * (r - hr)) / (2 * hr * r) \
           - (gr(r, th + hth) - gr(r, th - hth)) / (2 * hth * r)
    scale = a * r ** (m - 1.0)
    return np.max(np.abs(eq)) / scale, np.max(np.abs(comp)) / (gamma0 * (a * r ** (m - 1) / tau0) ** n)

eqr, compr = residuals()
print(f"   equilibrium residual (rel): {eqr:.3e};  compatibility residual (rel): {compr:.3e}")

# (c) constitutive: gamma parallel to tau and |gamma| = gamma0 (tau/tau0)^n; w odd, faces ok
th = np.linspace(0, np.pi, 20001)
tr = a * 0.03 ** (m - 1) * gp(th); tth = -a * m * 0.03 ** (m - 1) * g(th)
gr = gamma0 / tau0 ** n * a ** n * 0.03 ** ((m - 1) * n) * rho(th) ** (n - 1) * gp(th)
gth = -gamma0 / tau0 ** n * a ** n * 0.03 ** ((m - 1) * n) * m * rho(th) ** (n - 1) * g(th)
cross = tr * gth - tth * gr
print(f"   max |tau x gamma| / |tau||gamma| = {np.max(np.abs(cross)) / np.max(np.hypot(tr, tth) * np.hypot(gr, gth)):.3e}")
print(f"   face traction tau_thz(pi) = {tau_th(0.03, np.pi):.3e} MPa;  w(0.03, pi) = {w_disp(0.03, np.pi):.3e} mm (upper face +z: positive)")
