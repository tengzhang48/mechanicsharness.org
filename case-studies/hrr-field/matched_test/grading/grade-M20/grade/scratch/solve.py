import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import CubicSpline

# --- parameters ---
n = 3.5
m = 1.0/n
beta = 1.0/(n+1.0)
tau0 = 80.0      # MPa
gamma0 = 0.002
Jtarget = 0.120  # N/mm

# --- angular ODE:  (R f')' + D beta R f = 0,  R = (beta^2 f^2 + f'^2)^((m-1)/2)
# with D = 1 + (beta-1)/n = n/(n+1)
# state y = [f, u=f']
D = 1.0 + (beta-1.0)/n

def rhs(theta, y):
    f, u = y
    Q2 = beta*beta*f*f + u*u
    denom = beta*beta*f*f + m*u*u
    num = beta*f*((m-1.0)*beta*u*u + D*(beta*beta*f*f + u*u))
    up = -num/denom
    return [u, up]

sol = solve_ivp(rhs, [0.0, np.pi], [0.0, 1.0], rtol=1e-12, atol=1e-14,
                method='DOP853', dense_output=True, max_step=0.01)

print("status", sol.status, sol.message)
fpi = sol.y[0, -1]
upi = sol.y[1, -1]
print("f(pi) =", fpi)
print("f'(pi) =", upi, " (should be ~0)")

# evaluate f, f' on [-pi, pi] using odd symmetry: f(-t) = -f(t)
def sgn(th):
    th = np.asarray(th, dtype=float)
    return np.where(th < 0.0, -1.0, 1.0)

def fval(th):
    th = np.asarray(th, dtype=float)
    return sgn(th)*sol.sol(np.abs(th))[0]

def uval(th):
    th = np.asarray(th, dtype=float)
    return sgn(th)*sol.sol(np.abs(th))[1]

def Qval(th):
    th = np.asarray(th, dtype=float)
    f = fval(th); u = uval(th)
    return np.sqrt(beta*beta*f*f + u*u)

# --- J integral at A=1 on a circle radius r ---
# gamma_r = A beta r^{beta-1} f ; gamma_theta = A r^{beta-1} f'
# |gamma| = A r^{beta-1} Q
# tau = tau0 * gamma0^{-1/n} * (A r^{beta-1} Q)^{1/n}
# tau_r = tau * beta f / Q
# w_,x = gamma_r cos - gamma_theta sin
# W = n/(n+1) * tau * |gamma|

def J_of(A, r):
    th = np.linspace(-np.pi, np.pi, 200001)
    f = fval(th); u = uval(th)
    Q = Qval(th)
    gmag = A*r**(beta-1.0)*Q
    tau = tau0 * gamma0**(-1.0/n) * gmag**(1.0/n)
    gr = A*beta*r**(beta-1.0)*f
    gt = A*r**(beta-1.0)*u
    W = (n/(n+1.0))*tau*gmag
    taux = tau*gr/gmag
    w_x = gr*np.cos(th) - gt*np.sin(th)
    taur = tau*gr/gmag
    integrand = W*np.cos(th) - taur*w_x
    return np.trapezoid(integrand, th)*r

J1 = J_of(1.0, 1.0)
J2 = J_of(1.0, 2.0)
J3 = J_of(1.0, 0.02)
print("J(A=1,r=1) =", J1)
print("J(A=1,r=2) =", J2, " ratio", J2/J1)
print("J(A=1,r=0.02)=", J3, " ratio", J3/J1)

# J = alpha * A^{1+1/n}
alpha = J1
A = (Jtarget/alpha)**(1.0/(1.0+1.0/n))
print("A =", A)

# check J with A
print("J(A,r=1) =", J_of(A,1.0), " target", Jtarget)

# --- requested quantities ---
r1, r2 = 0.0200, 0.0800
Dw1 = 2.0*A*r1**beta*fpi
Dw2 = 2.0*A*r2**beta*fpi
print("[1] slip at 0.02:", Dw1, " at 0.08:", Dw2)

# stresses at r=0.02, theta=60deg
rq = 0.0200
thq = np.deg2rad(60.0)
f = float(fval(thq)); u = float(uval(thq)); Q = float(Qval(thq))
gmag = A*rq**(beta-1.0)*Q
tau = tau0*gamma0**(-1.0/n)*gmag**(1.0/n)
trz = tau*(beta*f)/Q
ttz = tau*(u)/Q
print("[2] tau_rz, tau_theta z at 0.02,60deg:", trz, ttz, " |gamma|=", gmag)

# radii where |gamma|=0.01 for theta=0,90,180
gtarget = 0.0100
for deg in [0.0, 90.0, 180.0]:
    th = np.deg2rad(deg)
    Q = float(Qval(th))
    r = (A*Q/gtarget)**(1.0/(1.0-beta))
    print("[3] |gamma|=0.01 at theta=%g deg -> r =" % deg, r)

# boundary residual: tau_theta z at faces
for th in [np.pi-1e-9, -np.pi+1e-9]:
    f=float(fval(th));u=float(uval(th));Q=float(Qval(th))
    gmag=A*1.0**(beta-1)*Q
    tau=tau0*gamma0**(-1.0/n)*gmag**(1.0/n)
    print("face theta=",th,"tau_theta z=",tau*u/Q)
