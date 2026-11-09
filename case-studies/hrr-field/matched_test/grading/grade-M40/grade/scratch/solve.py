"""Mode III (anti-plane) power-law crack-tip field, amplitude fixed by J.

Model: small strain, quasistatic, stationary slit along -x, material -pi<theta<pi.
w = w(x,y); gamma = grad w (engineering shear); law  gamma = (g0/t0)(tau/t0)^(n-1) tau,
inverted  tau = C |gamma|^(k-1) gamma  with k = 1/n, C = t0*g0^(-k).
Equilibrium div tau = 0.

Separable leading field:  w = A r^a f(theta),  a = 1/(n+1)  (finite nonzero J fixes it),
f(0)=0 (antisymmetry), f'(0)=1 (normalization).
ODE on [0, pi]:  (h^(k-1) f')' + (n/(n+1)) a f h^(k-1) = 0,  h = sqrt(a^2 f^2 + f'^2).
Traction-free faces theta=+-pi need tau_th=0 -> f'(pi)=0 (CHECKED, not imposed).

Units: mm, N, MPa.
"""
import numpy as np
from scipy.integrate import solve_ivp

# ---- material / loading ------------------------------------------------------
n    = 3.5
t0   = 80.0        # MPa
g0   = 0.002
J    = 0.120       # N/mm
k    = 1.0 / n
a    = 1.0 / (n + 1.0)          # displacement exponent; gamma ~ r^(a-1), tau ~ r^((a-1)/n)
c    = n / (n + 1.0)            # W = (n/(n+1)) tau*gamma (pure power law)
Cmat = t0 * g0 ** (-k)          # tau = C |gamma|^(k-1) gamma  [MPa]

# ---- angular ODE: state y = (f, fp) ------------------------------------------
def rhs(theta, y):
    f, fp = y
    h2 = a * a * f * f + fp * fp
    num = -c * a * f - (k - 1.0) * a * a * f * fp * fp / h2
    den = 1.0 + (k - 1.0) * fp * fp / h2
    return (fp, num / den)

sol = solve_ivp(rhs, (0.0, np.pi), (0.0, 1.0), rtol=1e-12, atol=1e-14,
                dense_output=True, max_step=0.01)

def f_of(t):
    t = np.asarray(t, float)
    return np.sign(t) * sol.sol(np.abs(t))[0]

def fp_of(t):
    t = np.asarray(t, float)
    return sol.sol(np.abs(t))[1]          # f' is even in theta

th = np.linspace(0.0, np.pi, 20001)
f, fp = sol.sol(th)
h  = np.sqrt(a * a * f * f + fp * fp)
face_residual = abs(fp[-1]) / (a * abs(f[-1]))

# ---- angular J integral ------------------------------------------------------
# J = C A^(k+1) I_th;  all r-factors cancel because (a-1)(k+1) = -1.
TH = np.linspace(-np.pi, np.pi, 400001)
F, FP = f_of(TH), fp_of(TH)
HH = np.sqrt(a * a * F * F + FP * FP)
integrand = (c * HH**(k + 1.0) * np.cos(TH)
             - HH**(k - 1.0) * a * F * (a * F * np.cos(TH) - FP * np.sin(TH)))
I_th = np.trapezoid(integrand, TH)

A = (J / (Cmat * I_th)) ** (1.0 / (k + 1.0))       # w = A r^a f

# ---- field evaluators --------------------------------------------------------
def w_field(x, y):
    r = np.hypot(x, y); t = np.arctan2(y, x)
    return A * r**a * f_of(t)

def gamma_polar(r, t):
    g = A * r**(a - 1.0)
    return g * a * f_of(t), g * fp_of(t)

def tau_polar(r, t):
    F_, FP_ = f_of(t), fp_of(t)
    HH_ = np.sqrt(a * a * F_ * F_ + FP_ * FP_)
    amp = Cmat * A**k * r ** ((a - 1.0) * k) * HH_**(k - 1.0)
    return amp * a * F_, amp * FP_

def tau_cart(x, y):
    r = np.maximum(np.hypot(x, y), 1e-15); t = np.arctan2(y, x)
    tr, tt = tau_polar(r, t)
    ct, st = np.cos(t), np.sin(t)
    return tr * ct - tt * st, tr * st + tt * ct

def W_den(x, y):
    gr, gt = gamma_polar(np.hypot(x, y), np.arctan2(y, x))
    g = np.hypot(gr, gt)
    return c * (Cmat * g**k) * g

def dw_dx(x, y):
    d = 1e-7
    return (w_field(x + d, y) - w_field(x - d, y)) / (2 * d)

# ---- requested quantities ----------------------------------------------------
print("=== solver diagnostics ===")
print(f"converged: {sol.success}; face traction residual |tau_th(pi)|/(a f(pi)) = {face_residual:.3e}")
print(f"f(pi) = {f[-1]:.10f},  I_th = {I_th:.10f},  A = {A:.10e} mm^(7/9)")

print("\n=== 1. relative face slip Delta w = w(r,+pi) - w(r,-pi) ===")
for r in (0.0200, 0.0800):
    print(f"  Delta w({r:.4f}) = {2*A*r**a*f[-1]:.6e} mm")

print("\n=== 2. stresses at r = 0.0200 mm, theta = 60 deg ===")
r0, th0 = 0.0200, np.deg2rad(60.0)
tr, tt = tau_polar(r0, th0)
gr, gt = gamma_polar(r0, th0)
print(f"  tau_rz = {float(tr):.6f} MPa,  tau_thz = {float(tt):.6f} MPa  (|gamma| = {float(np.hypot(gr,gt)):.4e})")

print("\n=== 3. radii where |gamma| = 0.0100 ===")
for name, t in (("0", 0.0), ("90 deg", np.pi / 2), ("180 deg", np.pi)):
    hh = float(np.sqrt(a*a*f_of(t)**2 + fp_of(t)**2))
    print(f"  theta = {name}: r = {(A*hh/0.0100)**(1.0/(1.0-a)):.6e} mm   [h(theta) = {hh:.6f}]")

# ---- CHECK A: div tau = 0 by Cartesian finite differences --------------------
print("\n=== check A: Cartesian equilibrium residual ===")
d = 1e-6
for (x0, y0) in ((0.02*np.cos(np.pi/3), 0.02*np.sin(np.pi/3)),
                 (0.05, 0.05), (-0.03, 0.04), (0.01, -0.008), (-0.06, 0.005)):
    dtx = (tau_cart(x0+d, y0)[0] - tau_cart(x0-d, y0)[0])/(2*d)
    dty = (tau_cart(x0, y0+d)[1] - tau_cart(x0, y0-d)[1])/(2*d)
    tx, ty = tau_cart(x0, y0)
    rel = abs(dtx + dty)/(abs(tx)/abs(x0) + abs(ty)/abs(y0) + 1e-300)
    print(f"  (x,y)=({x0:+.3f},{y0:+.3f}): div tau = {dtx+dty:+.3e}, relative {rel:.2e}")

# ---- CHECK B: J on closed rectangular contours (Cartesian derivatives) -------
print("\n=== check B: J on closed rectangles (target 0.120 N/mm) ===")
def J_rect(cx, cy, half, m=2001):
    """CCW rectangle; outward normal (dy,-dx) of each segment."""
    xs, ys = np.linspace(cx-half, cx+half, m), np.linspace(cy-half, cy+half, m)
    segs = [(xs[0], cy-half, xs[-1], cy-half),   # bottom -> +x
            (cx+half, ys[0], cx+half, ys[-1]),   # right  -> +y
            (cx+half, cy+half, cx-half, cy+half),# top    -> -x
            (cx-half, cy+half, cx-half, cy-half)]# left   -> -y
    total = 0.0
    for x1, y1, x2, y2 in segs:
        u = np.linspace(0.0, 1.0, 2001)
        xx, yy = x1 + (x2-x1)*u, y1 + (y2-y1)*u
        L = np.hypot(x2-x1, y2-y1)
        Nx, Ny = (y2-y1)/L, -(x2-x1)/L
        txz, tyz = tau_cart(xx, yy)
        tz = txz*Nx + tyz*Ny
        integ = W_den(xx, yy)*Nx - tz*dw_dx(xx, yy)
        total += np.trapezoid(integ, u)*L
    return total

for half in (0.05, 0.12):
    print(f"  square half-width {half:.2f} mm centered at tip: J = {J_rect(0.0, 0.0, half):.6f} N/mm")

# ---- CHECK C: n=1 elastic limit vs Westergaard -------------------------------
print("\n=== check C: n=1 limit ===")
a_, k_, c_ = 0.5, 1.0, 0.5
def rhs1(theta, y):
    return (y[1], -0.25*y[0])
s1 = solve_ivp(rhs1, (0, np.pi), (0.0, 1.0), rtol=1e-12, atol=1e-14, dense_output=True)
t1 = np.linspace(0, np.pi, 20001); f1, fp1 = s1.sol(t1)
print(f"  max |f - 2 sin(theta/2)| = {np.max(np.abs(f1 - 2*np.sin(t1/2))):.2e}")
THf = np.linspace(-np.pi, np.pi, 400001)
Ff  = np.sign(THf)*s1.sol(np.abs(THf))[0]; FPf = s1.sol(np.abs(THf))[1]
HHf = np.sqrt(0.25*Ff**2 + FPf**2)
integf = 0.5*HHf**2*np.cos(THf) - 0.5*Ff*(0.5*Ff*np.cos(THf) - FPf*np.sin(THf))
I1 = np.trapezoid(integf, THf)
G_ = t0/g0
A1 = (J/(G_*I1))**0.5
# elastic: w = (K/G) sqrt(2r/pi) sin(th/2) = A1 r^1/2 * 2 sin(th/2)  ->  K = A1 G sqrt(2 pi)
K_implied = A1*G_*np.sqrt(2*np.pi)
K_exact = np.sqrt(2*G_*J)
print(f"  I_th(n=1) = {I1:.10f}  (closed form: pi = {np.pi:.10f})")
print(f"  K_III from field = {K_implied:.8f}, exact sqrt(2 G J) = {K_exact:.8f}, rel {abs(K_implied/K_exact-1):.2e}")

# ---- CHECK D: w_,x consistency: FD vs analytic gamma_x ------------------------
print("\n=== check D: dw/dx FD vs analytic gamma_x ===")
for (x0, y0) in ((0.02*np.cos(np.pi/3), 0.02*np.sin(np.pi/3)), (-0.03, 0.04)):
    fd = float(dw_dx(x0, y0))
    t = np.arctan2(y0, x0); r = np.hypot(x0, y0)
    gr, gt = gamma_polar(r, t)
    an = float(gr*np.cos(t) - gt*np.sin(t))
    print(f"  (x,y)=({x0:+.3f},{y0:+.3f}): FD {fd:+.8f}, analytic {an:+.8f}, rel {abs(fd/an-1):.2e}")
