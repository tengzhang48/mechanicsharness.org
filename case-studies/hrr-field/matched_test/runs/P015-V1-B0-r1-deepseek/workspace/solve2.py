import numpy as np
from scipy.integrate import solve_ivp

tau0 = 80.0; gamma0 = 0.002; n = 3.5; J = 0.120
p = 1.0/(n+1.0); s = 1.0/(n+1.0)
C = tau0 * gamma0**(-1.0/n)

def rhs(theta, y):
    f, fp = y
    g2 = (s*f)**2 + fp**2
    g = np.sqrt(max(g2, 1e-300))
    k = 1.0/n - 1.0
    A22 = g**k + k * g**(k-2.0) * fp**2
    b = (1.0-p) * g**k * s * f + k * g**(k-2.0) * s**2 * f * fp**2
    return [fp, -b/A22]

sol = solve_ivp(rhs, [0.0, np.pi], [0.0, 1.0], rtol=1e-13, atol=1e-15, dense_output=True, method='DOP853')

def fth(th): return sol.sol(th)[0]
def fph(th): return sol.sol(th)[1]
def gth(th):
    f, fp = sol.sol(th); return np.sqrt((s*f)**2 + fp**2)

# --- J normalization integral (A=1), integrand coefficient of r^0 ---
def integrand(th):
    f, fp = fth(th), fph(th)
    g = gth(th)
    W1 = n/(n+1.0) * C * g**((n+1.0)/n)            # W = W1 r^-1
    Tr = C * g**(1.0/n-1.0) * (s*f)                # tau_r = Tr r^-p
    Wx1 = s*f*np.cos(th) - fp*np.sin(th)           # w_,x = Wx1 r^(s-1)
    return W1*np.cos(th) - Tr*Wx1

th = np.linspace(-np.pi, np.pi, 200001)
IJ = np.trapezoid(integrand(th), th)
print("I_J =", IJ)
A = (J/IJ)**(n/(n+1.0))
print("A =", A)

# independent check: build full physical field and integrate J on circles
def gamma_vec(r, th):
    return np.array([A*r**(s-1)*s*fth(th), A*r**(s-1)*fph(th)])  # (gamma_r, gamma_th)
def tau_vec(r, th):
    gr, gt = gamma_vec(r, th)
    g = np.hypot(gr, gt)
    return C*g**(1.0/n-1.0)*np.array([gr, gt])
def Jnum(r):
    a = np.linspace(-np.pi, np.pi, 200001)
    out = np.zeros_like(a)
    for i, t in enumerate(a):
        gr, gt = gamma_vec(r, t)
        tr, tt = tau_vec(r, t)
        W = n/(n+1.0)*(tr*gr+tt*gt)
        g = np.hypot(gr, gt)
        wx = gr*np.cos(t) - gt*np.sin(t)
        out[i] = (W*np.cos(t) - tr*wx)*r
    return np.trapezoid(out, a)
for r in [0.02, 0.05, 0.08, 0.2]:
    print(f"J numeric at r={r}: {Jnum(r):.8f}")

print()
# (1) face slip
for r in [0.0200, 0.0800]:
    w_face = A*r**s*fth(np.pi)
    print(f"dw(r={r}) = {2*w_face:.8f} mm")

# (2) stresses at r=0.02, theta=60deg
deg = np.pi/180
r0, th0 = 0.0200, 60*deg
gr, gt = gamma_vec(r0, th0)
tr, tt = tau_vec(r0, th0)
print(f"tau_rz(r=0.02,60) = {tr:.8f} MPa ; tau_thz = {tt:.8f} MPa")
print(f"|gamma| there = {np.hypot(gr,gt):.8f}")

# (3) radii where |gamma|=0.01
for name, t in [("0deg",0.0), ("90deg",np.pi/2), ("180deg",np.pi)]:
    g = gth(t)
    r = (0.01/(A*g))**(1.0/(s-1.0))
    print(f"|gamma|=0.01 at theta={name}: r={r:.8f} mm")
