"""
Mode III (anti-plane) stationary slit in a pure power-law shear material.

Constitutive law (given):  gamma = (gamma0/tau0) (tau/tau0)^(n-1) tau
  <=>  tau = tau0 (gamma/gamma0)^(1/n) * unit(gamma)     (nonlinear "p-harmonic" material)
Equilibrium: div tau = 0,  gamma = grad w.

Governing near-tip model: separable (HRR-type) field fixed by a finite, nonzero J:
  J ~ r^(1 - (alpha_w-1)(1+1/n))  =>  alpha_w = 1/(n+1),  w ~ r^(1/(n+1)),
  tau ~ r^(-1/(n+1)) (stress exponent p = 1/(n+1), the mode III HRR exponent).

Closed form (derived here, verified below): use the stress-function angle psi with
t = tan(psi) in [0, inf) for the upper half plane 0 <= theta <= pi:
  theta(t)  = atan(t) + atan(n t)                    (theta: 0 -> pi)
  S_r(t)    = a t / sqrt(1+n^2 t^2),  S_th(t) = a / sqrt(1+n^2 t^2),  a = n/(n+1)
  T(t)      = a sqrt((1+t^2)/(1+n^2 t^2))            (T^2 = S_r^2 + S_th^2)
  tau_rz = A r^(a-1) S_r,   tau_thz = A r^(a-1) S_th          (a-1 = -1/(n+1))
  |gamma| = c A^n r^(-n/(n+1)) T^n,   c = gamma0/tau0^n
  gamma_r = |gamma| S_r/T,   gamma_th = |gamma| S_th/T
  w = D r^(1/(n+1)) F(theta),  D = c A^n,  F(theta) = int_0^theta T^n cos(psi') dtheta'
n = 1 limit: theta = 2 psi, g = -cos(theta/2) etc => exact LEFM field.

Amplitude A fixed by J = 0.120 N/mm (J scales as A^(n+1)).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

# ---------------- parameters
tau0, gamma0, n = 80.0, 0.002, 3.5
J_target = 0.120                      # N/mm
cc = gamma0 / tau0**n                  # strain per (MPa)^n
alpha_w = 1.0 / (n + 1.0)              # displacement exponent of w
a_exp = n / (n + 1.0)                  # stress-function exponent: tau ~ r^(a_exp-1)

# ---------------- closed-form angular functions, tabulated vs theta
t_grid = np.concatenate([np.linspace(0.0, 1.0, 4001),
                         np.geomspace(1.0, 1e9, 12001)[1:]])
th_t = np.arctan(t_grid) + np.arctan(n * t_grid)          # theta(t)
a = n / (n + 1.0)
T_t = a * ((1 + t_grid**2) / (1 + n**2 * t_grid**2))**(1.0/(2*(n+1.0)))
S_r_t = T_t * t_grid / np.sqrt(1 + t_grid**2)             # = g' = T sin(psi)
S_th_t = T_t / np.sqrt(1 + t_grid**2)                     # = -a g = T cos(psi)
cosp = 1.0 / np.sqrt(1 + t_grid**2)                       # cos(psi)
dth_dt = 1.0/(1 + t_grid**2) + n/(1 + n**2 * t_grid**2)   # dtheta/dt
F_t = np.concatenate([[0.0], np.cumsum(
    0.5*(T_t[:-1]**n * cosp[:-1] * dth_dt[:-1] + T_t[1:]**n * cosp[1:] * dth_dt[1:]) * np.diff(t_grid))])

# enforce monotone theta for splines
order = np.argsort(th_t)
S_r_s  = CubicSpline(th_t[order], S_r_t[order])
S_th_s = CubicSpline(th_t[order], S_th_t[order])
T_s    = CubicSpline(th_t[order], T_t[order])
F_s    = CubicSpline(th_t[order], F_t[order])
th_max = th_t[-1]                                         # = pi (to machine precision)

def t_of_theta(th):
    f = lambda t: np.arctan(t) + np.arctan(n*t) - th
    hi = 1.0
    while f(hi) < 0: hi *= 2
    return brentq(f, 0.0, hi, xtol=1e-14, rtol=1e-15)

def field(r, th, A):
    """returns w, tau_rz, tau_thz, |gamma|, gamma_r, gamma_th (th in [0,pi])"""
    Sr, St, T, F = S_r_s(th), S_th_s(th), T_s(th), F_s(th)
    D = cc * A**n
    w = D * r**alpha_w * F
    tau_r, tau_th = A * r**(a_exp-1) * Sr, A * r**(a_exp-1) * St
    gm = D * r**(alpha_w-1) * T**n
    return w, tau_r, tau_th, gm, gm*Sr/T, gm*St/T

# ---------------- J of the field on a circular arc (Cartesian components)
def J_circle(R, A, k=200001):
    th = np.linspace(-np.pi, np.pi, k)
    th_u = np.abs(th)                     # symmetric angular functions
    Sr, St, T = S_r_s(th_u), S_th_s(th_u), T_s(th_u)
    D = cc * A**n
    tau_r = A * R**(a_exp-1) * Sr
    tau_th = A * R**(a_exp-1) * St * np.sign(th)          # odd in theta
    gmag = D * R**(alpha_w-1) * T**n
    gam_r = gmag * Sr / T
    gam_th = gmag * St / T * np.sign(th)                  # odd
    W = (n/(n+1.0)) * (A * R**(a_exp-1) * T) * gmag       # = n/(n+1) * tau.mag * gamma.mag
    gam_x = gam_r*np.cos(th) - gam_th*np.sin(th)      # Cartesian component of dw/dx
    integrand = W*np.cos(th) - tau_r*gam_x
    return np.trapezoid(integrand, th) * R

# amplitude from J  (J = A^(n+1) * J(A=1))
J1 = J_circle(0.05, A=1.0)
A_amp = (J_target / J1)**(1.0/(n+1.0))
D_amp = cc * A_amp**n

print("=== amplitude ===")
print(f"J(A=1) at R=0.05      = {J1:.8e} N/mm (units with A=1)")
print(f"A (tau = A r^(-1/(n+1)) S)  = {A_amp:.8f} MPa mm^(1/(n+1))")
print(f"D = c A^n (w = D r^(1/(n+1)) F) = {D_amp:.8e} mm^(n/(n+1))")

# J path-independence with final amplitude
print("\n=== verification 1: J on several circles ===")
for R in (0.005, 0.02, 0.08, 0.20, 0.30):
    print(f"  R = {R:6.3f} mm:  J = {J_circle(R, A_amp):.8f} N/mm")

# ---------------- answers
print("\n=== answers ===")
for r in (0.0200, 0.0800):
    w_up, *_ = field(r, np.pi, A_amp)
    print(f"1) dw(r={r:.4f}) = {2*w_up:.6e} mm")

r2, th2 = 0.0200, np.deg2rad(60.0)
_, tr, tt, *_ = field(r2, th2, A_amp)
print(f"2) r=0.02, th=60deg: tau_rz = {tr:.4f} MPa, tau_thz = {tt:.4f} MPa")

gmag_tgt = 0.0100
for thd in (0.0, 90.0, 180.0):
    th = np.deg2rad(thd)
    T = float(T_s(th)) if thd < 180.0 else a*n**(-1.0/(n+1.0))  # T(pi) = a n^(-1/(n+1))
    r_g = (D_amp * T**n / gmag_tgt)**(1.0/(1.0-alpha_w))
    print(f"3) |gamma|=0.01 at theta={thd:5.1f} deg: r = {r_g:.6f} mm  (T={T:.6f})")

# ---------------- verification 2: n=1 regression against exact LEFM
print("\n=== verification 2: n=1 regression (exact LEFM) ===")
def angular_splines(nv):
    tg = np.concatenate([np.linspace(0, 1, 4001), np.geomspace(1, 1e6, 8001)[1:]])
    th = np.arctan(tg) + np.arctan(nv*tg)
    av = nv/(nv+1)
    T = av*((1+tg**2)/(1+nv**2*tg**2))**(1.0/(2*(nv+1)))
    Sr = T*tg/np.sqrt(1+tg**2); St = T/np.sqrt(1+tg**2)
    cp = 1/np.sqrt(1+tg**2); dth = 1/(1+tg**2)+nv/(1+nv**2*tg**2)
    Ft = np.concatenate([[0], np.cumsum(0.5*(T[:-1]**nv*cp[:-1]*dth[:-1]+T[1:]**nv*cp[1:]*dth[1:])*np.diff(tg))])
    o = np.argsort(th)
    return CubicSpline(th[o], Sr[o]), CubicSpline(th[o], St[o]), CubicSpline(th[o], T[o]), CubicSpline(th[o], Ft[o])

Sr1, St1, T1, F1 = angular_splines(1.0)
thc = np.linspace(0, np.pi, 2001)
print("  max|F(theta) - sin(theta/2)| =", np.max(np.abs(F1(thc) - np.sin(thc/2))))
print("  max|S_r - 0.5 sin(th/2)| =", np.max(np.abs(Sr1(thc) - 0.5*np.sin(thc/2))),
      "  max|S_th - 0.5 cos(th/2)| =", np.max(np.abs(St1(thc) - 0.5*np.cos(thc/2))))
# LEFM amplitude check: with n=1, c=1/G, w = (A/G) r^1/2 sin(th/2) = (2K/G)sqrt(r/2pi) sin(th/2)
# => K = A sqrt(pi/2);  J = K^2/(2G) = A^2 pi/(4G)  =>  I_J(=J/(c A^2)) must equal pi/4
def J_circle_n1(R, A):
    th = np.linspace(-np.pi, np.pi, 200001); tu = np.abs(th)
    Sr, St, T = Sr1(tu), St1(tu), T1(tu)
    G = tau0/gamma0
    tau_r = A*R**-0.5*Sr; tau_th = A*R**-0.5*St*np.sign(th)
    gm = (1.0/G)*A*1.0*R**-0.5*T
    gr = gm*Sr/T; gth = gm*St/T*np.sign(th)
    W = 0.5*(A*R**-0.5*T)*gm
    gx = gr*np.cos(th)-gth*np.sin(th)
    return np.trapezoid((W*np.cos(th)-tau_r*gx)*R, th)
J1n1 = J_circle_n1(0.05, 1.0)
print("  I_J(n=1) =", J1n1*G if False else J1n1/(gamma0/tau0), "  pi/4 =", np.pi/4)

# ---------------- verification 3: ODE cross-check of the angular functions
print("\n=== verification 3: nonlinear-oscillator ODE vs closed form ===")
# stream-function form: g(0)=-1, g'(0)=0,  (T^(n-1) g')' + k T^(n-1) g = 0,
#   T^2 = g'^2 + a^2 g^2,  k = n/(n+1)^2;  expect g(pi)=0 and g' matching closed form.
k_ode = n/(n+1.0)**2
def rhs(th, y):
    g, gp = y
    T2 = gp**2 + a**2*g**2
    gpp = -g*(k_ode*T2 + (n-1)*a**2*gp**2)/(T2 + (n-1)*gp**2)
    return [gp, gpp]
sol = solve_ivp(rhs, [0, np.pi], [-1.0, 0.0], rtol=1e-12, atol=1e-14,
                dense_output=True, max_step=np.pi/2000)
ths = np.linspace(0, np.pi, 41)
g_num = sol.sol(ths)[0]
def g_closed(th):
    t = t_of_theta(th)
    return -(((1+t**2)/(1+n**2*t**2))**(1.0/(2*(n+1))))/np.sqrt(1+t**2)
def gp_closed(th):
    t = t_of_theta(th)
    return a*t*((1+t**2)/(1+n**2*t**2))**(1.0/(2*(n+1)))/np.sqrt(1+t**2)
err_g = max(abs(g_num[i]-g_closed(th)) for i, th in enumerate(ths))
print(f"  g(pi) from ODE = {sol.y[0,-1]:.3e}  (should be 0: traction-free + J consistency)")
print(f"  max|g_ODE - g_closed| = {err_g:.3e}")
print(f"  g'(pi) ODE = {sol.y[1,-1]:.6f}, closed = {gp_closed(np.pi):.6f}  (= T(pi) = a n^(-1/(n+1)) = {a*n**(-1/(n+1)):.6f})")

# ---------------- verification 4: Cartesian FD equilibrium + constitutive + faces
print("\n=== verification 4: PDE/BC residuals of the full field ===")
def tau_cart(x, y):
    r = np.hypot(x, y); th = np.arctan2(y, x); tu = np.abs(th)
    Sr, St = S_r_s(tu), S_th_s(tu)
    tr = A_amp*r**(a_exp-1)*Sr; tt = A_amp*r**(a_exp-1)*St*np.sign(th)
    return tr*np.cos(th)-tt*np.sin(th), tr*np.sin(th)+tt*np.cos(th)
h = 1e-5
res_eq, mag_eq = [], []
for rr in (0.02, 0.1):
    for thd in (30.0, 90.0, 150.0):
        x, y = rr*np.cos(np.deg2rad(thd)), rr*np.sin(np.deg2rad(thd))
        dtxx = (tau_cart(x+h, y)[0]-tau_cart(x-h, y)[0])/(2*h)
        dtyy = (tau_cart(x, y+h)[1]-tau_cart(x, y-h)[1])/(2*h)
        res_eq.append(abs(dtxx+dtyy)); mag_eq.append(abs(dtxx)+abs(dtyy))
print(f"  max FD equilibrium residual / scale = {max(r/m for r, m in zip(res_eq, mag_eq)):.3e}")
print(f"  tau_thz on faces (theta=+-pi): {A_amp*0.05**(a_exp-1)*float(S_th_s(np.pi)):.3e} MPa at r=0.05 (exact 0 by construction)")
# gamma is a gradient: check grad w by FD against analytic gamma at the SAME point
def w_of(x, y):
    r = np.hypot(x, y); th = np.arctan2(y, x); tu = np.abs(th)
    return D_amp*r**alpha_w*F_s(tu)*np.sign(th)
r_, th_ = 0.03, np.deg2rad(45.0)
_, tr_, tt_, gm_, gr_, gth_ = field(r_, th_, A_amp)
tau_mag = np.hypot(tr_, tt_)
print(f"  constitutive check: |tau| = {tau_mag:.6f} MPa vs tau0(|g|/g0)^(1/n) = {tau0*(gm_/gamma0)**(1/n):.6f} MPa")
x0, y0 = r_*np.cos(th_), r_*np.sin(th_)
gx_fd = (w_of(x0+h, y0)-w_of(x0-h, y0))/(2*h)
gy_fd = (w_of(x0, y0+h)-w_of(x0, y0-h))/(2*h)
gx_an, gy_an = gr_*np.cos(th_)-gth_*np.sin(th_), gr_*np.sin(th_)+gth_*np.cos(th_)
print(f"  grad-w check at (r=0.03,45deg): FD ({gx_fd:.6e},{gy_fd:.6e}) vs analytic ({gx_an:.6e},{gy_an:.6e})")
# spot-check the reported |gamma|=0.01 contour radii
print("  contour spot checks (should be 0.01):")
for thd, rr in ((0.0, 0.046894), (90.0, 0.025066), (180.0, 0.013398)):
    print(f"    theta={thd:5.1f} r={rr}: |gamma| = {field(rr, np.deg2rad(min(thd,179.999)), A_amp)[3]:.6f}")
