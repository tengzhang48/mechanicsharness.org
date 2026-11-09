"""Independent, first-principles verification for n=6.5 plane strain.

Reconstructs the DIMENSIONAL HRR field from the angular solution with the task's
material constants and an arbitrary J, then evaluates everything (J, H, C)
directly: J by full-circle quadrature with Cartesian displacement gradients
obtained by finite differences of u(x,y); C by a genuine double integral of W.
"""
import sys, numpy as np
sys.path.insert(0, "knowledge")
import hrr

n = 6.5
p = 1.0 / (n + 1.0)
s = hrr.singular_exponent(n)

# material / loading (any positive J works: field is self-similar in J)
alpha, sig0, E = 0.02, 400.0e6, 200.0e9
eps0 = sig0 / E
Jval = 5.0e-5          # J/m^2
R = 1.0e-4             # m, contour radius

sol = hrr.solve_angular(n, plane_strain=True, n_theta=16001)
I = hrr.integral_I(n, plane_strain=True, sol=sol)
Amp = sig0 * (Jval / (alpha * eps0 * sig0 * I)) ** p      # A: sigma ~ A r^-p

# --- angular functions extended to [-pi, pi] by mode-I symmetry --------------
th_h = sol["theta"]                       # 0..pi grid
fpar = {"f": sol["f"], "fp": sol["fp"], "fpp": sol["fpp"], "fppp": sol["fppp"]}
def ang(theta):
    """f and derivatives at theta in [-pi,pi] (even/odd symmetry)."""
    z = np.abs(theta)
    sgn = np.sign(theta)
    idx = np.searchsorted(th_h, z)
    idx = np.clip(idx, 1, len(th_h) - 1)
    t0, t1 = th_h[idx - 1], th_h[idx]
    w = (z - t0) / (t1 - t0)
    out = []
    for k, key in enumerate(["f", "fp", "fpp", "fppp"]):
        y = fpar[key][idx - 1] * (1 - w) + fpar[key][idx] * w
        out.append(y * (sgn ** k) if np.isscalar(theta) else y * sgn ** k)
    return out

def stresses(theta):
    """polar stress components (dimensional) and sigma_e at radius-scaled r=1"""
    f, fp, fpp, _ = ang(theta)
    sr = s * f + fpp
    st = s * (s - 1.0) * f
    srt = (1.0 - s) * fp
    se = np.sqrt(0.75 * (sr - st) ** 2 + 3.0 * srt ** 2)   # plane strain
    return sr, st, srt, se

def utilde(theta):
    """angular displacement functions (u_r = aeps (A/sig0)^n r^p ur_t etc.)"""
    f, fp, fpp, fppp = ang(theta)
    sr, st, srt, se = stresses(theta)
    rho = se ** (n - 1.0)
    c1 = 0.75 * (sr - st)                    # = (3/2) s_rr deviator
    ur = rho * c1 / p                        # /p = /(n(s-2)+1)
    # derivatives by spectral-free central differences on the theta functions
    h = 1e-6
    fp_, _, _, _ = None, None, None, None
    urp = (utilde_raw(theta + h)[0] - utilde_raw(theta - h)[0]) / (2 * h)
    ut = (urp - 3.0 * rho * srt) / (1.0 - p)
    return ur, ut

_ur_store = {}
def utilde_raw(theta):
    f, fp, fpp, fppp = ang(theta)
    sr, st, srt, se = stresses(theta)
    rho = se ** (n - 1.0)
    return (rho * 0.75 * (sr - st) / p,)

def u_cart(x, y):
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    ur_t, ut_t = utilde(theta)
    scale = alpha * eps0 * (Amp / sig0) ** n * r ** p
    ur, ut = scale * ur_t, scale * ut_t
    return ur * np.cos(theta) - ut * np.sin(theta), ur * np.sin(theta) + ut * np.cos(theta)

# --- J by direct full-circle quadrature with FD Cartesian gradients ----------
Nth = 400001
thg = np.linspace(-np.pi, np.pi, Nth)
rg = R * np.ones_like(thg)
xg, yg = rg * np.cos(thg), rg * np.sin(thg)
dx = 1e-9
uxp = (u_cart(xg + dx, yg)[0] - u_cart(xg - dx, yg)[0]) / (2 * dx)
uyp = (u_cart(xg + dx, yg)[1] - u_cart(xg - dx, yg)[1]) / (2 * dx)

sr, st, srt, se = stresses(thg)
sig_rr = Amp * rg ** (-p) * sr
sig_tt = Amp * rg ** (-p) * st
sig_rt = Amp * rg ** (-p) * srt
W = (n / (n + 1.0)) * alpha * eps0 * sig0 * (Amp * rg ** (-p) * se / sig0) ** (n + 1.0)
tx = sig_rr * np.cos(thg) - sig_rt * np.sin(thg)
ty = sig_rr * np.sin(thg) + sig_rt * np.cos(thg)
integrand = (W * np.cos(thg) - (tx * uxp + ty * uyp)) * rg
J_direct = np.trapezoid(integrand, thg)

print(f"J given            = {Jval:.8e}")
print(f"J direct (full 2pi)= {J_direct:.8e}   rel diff = {J_direct/Jval-1:.3e}")

# --- H from the reconstructed field ------------------------------------------
sr0, st0, srt0, se0 = stresses(np.array(0.0))
sig_th0 = Amp * R ** (-p) * float(st0)
H_direct = sig_th0 / (sig0 * (Jval / (alpha * eps0 * sig0 * R)) ** p)
print(f"H from field       = {H_direct:.8f}")

# --- C from a genuine double integral of W ------------------------------------
Nr, Nt = 2000, 16001
rr = np.linspace(R / Nr / 2, R - R / Nr / 2, Nr)      # midpoint radii
th2 = np.linspace(-np.pi, np.pi, Nt)
TH, RR = np.meshgrid(th2, rr, indexing="ij")
sr2, st2, srt2, se2 = stresses(TH)
sig_e_dim = Amp * RR ** (-p) * se2
W2 = (n / (n + 1.0)) * alpha * eps0 * sig0 * (sig_e_dim / sig0) ** (n + 1.0)
U = np.trapezoid(np.trapezoid(W2 * RR, th2, axis=0), rr)
print(f"U(R) direct        = {U:.8e} J/m")
print(f"J*R                = {Jval*R:.8e}")
print(f"C = U/(JR) direct  = {U/(Jval*R):.8f}")

# --- p from the field itself (log-log slope) ----------------------------------
r_probe = np.array([1e-5, 3e-5, 1e-4, 3e-4, 1e-3])
sig_yy0 = Amp * r_probe ** (-p) * float(st0)
slope = np.polyfit(np.log(r_probe), np.log(sig_yy0), 1)[0]
print(f"log-log slope of sigma_th(r,0) = {slope:.10f}  -> p = {-slope:.10f}")
