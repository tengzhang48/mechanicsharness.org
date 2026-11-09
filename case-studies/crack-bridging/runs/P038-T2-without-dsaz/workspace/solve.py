import numpy as np
from scipy.optimize import fsolve, brentq

# ---------------- data ----------------
b    = 0.025            # width, m
E1   = 110e9; h1 = 0.60e-3   # upper strip
E2   = 190e9; h2 = 0.90e-3   # lower strip
a0   = 20.0e-3          # starter cut length
sig  = 0.20e6           # sigma_b, Pa
dc   = 0.20e-3          # delta_c, m
G0   = 1.40             # J/m^2

# full-width flexural rigidities
EI1 = E1*b*h1**3/12.0
EI2 = E2*b*h2**3/12.0
R   = 1.0/EI1 + 1.0/EI2          # 1/(N m^2)
print("EI1=%.6e  EI2=%.6e  R=%.6f" % (EI1, EI2, R))

# ---------------- energy: G = Gamma0 + cohesive work ----------------
Gc  = sig*dc/2.0                 # connector work per area (linear softening) = 20 J/m^2
G   = G0 + Gc
M   = np.sqrt(2.0*b*G/R)
print("cohesive work=%.4f J/m2, total G=%.4f J/m2" % (Gc, G))
print("M = %.6f N.m" % M)

# ---------------- cohesive zone shape (relative opening) ----------------
# delta'''' - k delta = -b sig R ;  particular delta_p = dc
k    = b*sig*R/dc
beta = k**0.25
print("beta = %.4f 1/m ; 2pi/beta = %.3f mm" % (beta, 1e3*2*np.pi/beta))

# general: delta(x) = dc + A cosh(bx) + B sinh(bx) + C cos(bx) + D sin(bx)
# tip x=0: delta=dc, delta''=M R, delta'''=0
# root x=L: delta=0, delta'=0
def delta_profile(x, A, B, C, D, ell):
    ux = beta*x
    return dc + A*np.cosh(ux) + B*np.sinh(ux) + C*np.cos(ux) + D*np.sin(ux)

# From tip conditions: A + C = 0 ; (A - C)*beta^2 = M*R ; B = D
A = M*R/(2.0*beta**2)
C = -A

# root conditions give B (and the length ell):
#   delta(L)=0 : A(cosh u - cos u) + B(sinh u + sin u) = -dc
#   delta'(L)=0: A(sinh u + sin u) + B(cosh u + cos u) = 0
def resid(ell):
    u = beta*ell
    # solve B from the second (linear) eq, then check first
    B = -A*(np.sinh(u)+np.sin(u))/(np.cosh(u)+np.cos(u))
    return A*(np.cosh(u)-np.cos(u)) + B*(np.sinh(u)+np.sin(u)) + dc

# scan for a sign change
lams = np.linspace(1e-4, 0.06, 2000)
vals = np.array([resid(l) for l in lams])
sc = np.where(np.diff(np.sign(vals)))[0]
print("sign changes at ell ~", lams[sc])
ell = brentq(resid, lams[sc[0]], lams[sc[0]+1])
u = beta*ell
B = -A*(np.sinh(u)+np.sin(u))/(np.cosh(u)+np.cos(u))
D = B
print("L = ell = %.6f mm" % (1e3*ell))

def delta_at(x):
    return delta_profile(x, A, B, C, D, ell)

print("delta(0)       = %.6f um" % (1e6*delta_at(0)))
print("delta(L/2)     = %.6f um" % (1e6*delta_at(ell/2)))
print("delta(L)       = %.6f um" % (1e6*delta_at(ell)))

# check monotone & first crossing of dc
xs = np.linspace(0, ell, 400)
ds = delta_at(xs)
print("max delta in zone = %.4f um at x=%.4f mm" % (1e6*ds.max(), 1e3*xs[np.argmax(ds)]))

# ---------------- alternative fully-determined BVP (delta''(L)=0) ----------------
def resid2(vars):
    Al, Bl, ell2 = vars
    u2 = beta*ell2
    def delv(x, A, B):
        ux=beta*x
        return dc + A*np.cosh(ux)+B*np.sinh(ux)-A*np.cos(ux)+B*np.sin(ux)
    d0  = delv(0,Al,Bl)                       # want dc
    dL  = delv(ell2,Al,Bl)                     # want 0
    d1  = A*np.sinh(u2)+B*np.cosh(u2)+A*np.sin(u2)+B*np.cos(u2)   # want 0
    d2  = A*np.cosh(u2)+B*np.sinh(u2)+A*np.cos(u2)-B*np.sin(u2)   # want 0
    return [d0-dc, dL, d1, d2]
sol = fsolve(lambda v: resid2(v), [1e-4, -1e-4, 0.012], full_output=False)
Al, Bl, ell2 = sol
u2 = beta*ell2
M2 = 2*Al*beta**2/R
print("\n[var2 delta''(L)=0] A=%.4e B=%.4e ell=%.4f mm  M=%.4f N.m" % (Al,Bl,1e3*ell2,M2))
def del2(x):
    ux=beta*x
    return dc+Al*np.cosh(ux)+Bl*np.sinh(ux)-Al*np.cos(ux)+Bl*np.sin(ux)
print("[var2] delta(L/2)=%.4f um, delta(0)=%.4f um" % (1e6*del2(ell2/2), 1e6*del2(0)))