import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ---- parameters ----
tau0 = 80.0
gamma0 = 0.002
n = 3.5
s = 1.0/(n+1.0)
alpha = (1.0-n)/n
q = (s-1.0)/n
Jimp = 0.120

def rhs(th, y):
    f, fp = y
    num = -f*((1+q)*s*s*s*f*f + ((1+q)*s + alpha*s*s)*fp*fp)
    den = s*s*f*f + fp*fp/n
    return [fp, num/den]

# integrate over full range with odd symmetry: start theta=0, f(0)=0, f'(0)=1
sol = solve_ivp(rhs, [0, np.pi], [0.0, 1.0], rtol=1e-13, atol=1e-15,
                dense_output=True, max_step=0.005)
fpi = sol.y[0,-1]; fppi = sol.y[1,-1]
print("f(pi) =", fpi, " f'(pi) =", fppi)

def f_fp(th):
    scalar = np.isscalar(th) or np.ndim(th) == 0
    th = np.atleast_1d(np.asarray(th, float))
    out_f = np.empty_like(th); out_fp = np.empty_like(th)
    pos = th >= 0
    if np.any(pos):
        fm = sol.sol(th[pos]); out_f[pos]=fm[0]; out_fp[pos]=fm[1]
    neg = ~pos
    if np.any(neg):
        fm = sol.sol(-th[neg]); out_f[neg]=-fm[0]; out_fp[neg]=fm[1]
    if scalar:
        return float(out_f[0]), float(out_fp[0])
    return out_f, out_fp

def gfun(th):
    f, fp = f_fp(th)
    return np.sqrt(s*s*f*f + fp*fp), f, fp

# ---- build the angular integral I_J ----
# J = tau0 * gamma0**(-1/n) * A**((n+1)/n) * I_J
# integrand = [n/(n+1) g**((n+1)/n) cos th - s f g**((1-n)/n)(s f cos th - f' sin th)]
th = np.linspace(-np.pi, np.pi, 2000001)
g, f, fp = gfun(th)
integrand = (n/(n+1))*g**((n+1)/n)*np.cos(th) - s*f*g**((1-n)/n)*(s*f*np.cos(th) - fp*np.sin(th))
I_J = np.trapezoid(integrand, th)
print("I_J =", I_J)

pref = tau0*gamma0**(-1.0/n)
A = (Jimp/(pref*I_J))**(n/(n+1))
print("A =", A)

def gamma_mag(r, th):
    g,_,_ = gfun(th)
    return A * r**(s-1) * g

def tau_vec(r, th):
    g, f, fp = gfun(th)
    gm = A*r**(s-1)*g
    gr = A*r**(s-1)*s*f
    gt = A*r**(s-1)*fp
    tmag = tau0*(gm/gamma0)**(1.0/n)
    tr = tmag*gr/gm
    tt = tmag*gt/gm
    return tr, tt, tmag, gm

# ---- Q1: face slip ----
for r in (0.02, 0.08):
    dw = 2*A*r**s*fpi
    print(f"Q1 r={r}: Delta w = {dw:.8g} mm")

# ---- Q2: stresses at r=0.02, theta=60deg ----
r0 = 0.0200; th0 = np.deg2rad(60.0)
tr, tt, tmag, gm = tau_vec(r0, th0)
fv, fpv = f_fp(th0)
print(f"Q2 r=0.02 theta=60: tau_rz={tr:.8g} MPa, tau_thz={tt:.8g} MPa, |tau|={tmag:.8g}, |gamma|={gm:.8g}")
print(f"   f={fv:.8g} f'={fpv:.8g} gamma_r={A*r0**(s-1)*s*fv:.8g} gamma_th={A*r0**(s-1)*fpv:.8g} g={np.sqrt(s*s*fv*fv+fpv*fpv):.8g}")

# ---- independent check: equilibrium div(sigma)=0 in Cartesian by finite differences ----
def stress_xy(x, y):
    r = np.hypot(x, y); th = np.arctan2(y, x)
    g, f, fp = gfun(th)
    gm = A*r**(s-1)*g
    tr = tau0*(gm/gamma0)**(1.0/n) * (s*f/g)
    tt = tau0*(gm/gamma0)**(1.0/n) * (fp/g)
    c, sn = np.cos(th), np.sin(th)
    return tr*c - tt*sn, tr*sn + tt*c

worst = 0.0
for deg in range(-175, 176, 5):
    a = np.deg2rad(deg); x, y = np.cos(a), np.sin(a)
    h = 1e-6
    dsxz_dx = (stress_xy(x+h, y)[0] - stress_xy(x-h, y)[0])/(2*h)
    dsyz_dy = (stress_xy(x, y+h)[1] - stress_xy(x, y-h)[1])/(2*h)
    worst = max(worst, abs(dsxz_dx + dsyz_dy))
print("max |div sigma| on unit circle (Cartesian FD) =", worst)
print("face shear |sigma_yz| at theta=pi =", abs(stress_xy(np.cos(np.pi), np.sin(np.pi+1e-12))[1]))


# ---- Q3: radii where |gamma|=0.01 ----
gval = 0.0100
for deg in (0.0, 90.0, 180.0):
    thv = np.deg2rad(deg)
    g,_,_ = gfun(thv)
    r = (gval/(A*g))**(1.0/(s-1))
    print(f"Q3 theta={deg}: r(|gamma|=0.01) = {r:.8g} mm, g(theta)={g:.8g}")

# ---- verification: J on several circles using full Cartesian field ----
def w_fun(x, y):
    r = np.hypot(x,y); th = np.arctan2(y,x)
    f,_ = f_fp(th)
    return A*r**s*f

def J_circle(R, k=400001):
    a = np.linspace(-np.pi, np.pi, k)
    x, y = R*np.cos(a), R*np.sin(a)
    th = a
    # stresses in polar -> cartesian for the antiplane vector
    g, f, fp = gfun(th)
    gm = A*R**(s-1)*g
    gr = A*R**(s-1)*s*f
    gt = A*R**(s-1)*fp
    tmag = tau0*(gm/gamma0)**(1.0/n)
    tr = tmag*gr/gm; tt = tmag*gt/gm
    c, sn = np.cos(th), np.sin(th)
    sxz = tr*c - tt*sn
    syz = tr*sn + tt*c
    # W = n/(n+1) tau*gamma
    W = n/(n+1)*tmag*gm
    # w_,x numerically
    h = 1e-6*R
    dwx = (w_fun(x+h,y)-w_fun(x-h,y))/(2*h)
    Nx, Ny = c, sn
    tz = sxz*Nx + syz*Ny
    integrand = (W*Nx - tz*dwx)*R
    return np.trapezoid(integrand, a)

for R in (0.5, 1.0, 2.0):
    print(f"J (R={R}) = {J_circle(R):.10g}")
