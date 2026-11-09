import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

# ---------------- constants (SI) ----------------
b    = 25.0e-3          # strip width [m]
E1,h1= 110e9, 0.60e-3   # upper strip
E2,h2= 190e9, 0.90e-3   # lower strip
d1   = E1*h1**3/12      # bending stiffness per unit width [N m]
d2   = E2*h2**3/12
dbar = d1*d2/(d1+d2)    # series combination: 1/dbar = 1/d1 + 1/d2
G0   = 1.40             # intrinsic fracture energy [J/m^2]
sb   = 0.20e6           # bridging strength [Pa]
dc   = 0.20e-3          # connector critical opening [m]

print(f"d1 = {d1:.6f} N m, d2 = {d2:.6f} N m, dbar = {dbar:.8f} N m")

# ---------------- cohesive-zone beam model ----------------
# delta = dc*(1-u);  u'''' = lam^4 u  on 0<x<L
lam  = (sb/(dc*dbar))**0.25
# tip condition: G = (dbar*delta''(0))^2/(2*dbar) = G0  =>  u''(0) = -sqrt(2*G0/dbar)/dc
Kcur = np.sqrt(2*G0/dbar)/dc
k    = Kcur/lam**2        # = sqrt(2*G0/(dc*sb))
print(f"lam = {lam:.6f} 1/m  (1/lam = {1e3/lam:.5f} mm),  k = sqrt(2G0/(sb*dc)) = {k:.9f}")

A = (1-k)/2
C = (1+k)/2
# u(x) = A cosh(lx) + B sinh(lx) + C cos(lx) - B sin(lx)
# BCs at s = lam L:  u=0 and u'''=0

def beta_of(s):
    return -(A*np.cosh(s) + C*np.cos(s))/(np.sinh(s) - np.sin(s))

def g(s):   # shear condition u'''(L)=0 after eliminating B with u(L)=0
    be = beta_of(s)
    return A*np.sinh(s) + be*(np.cosh(s) + np.cos(s)) + C*np.sin(s)

s = brentq(g, 1.5, 1.8, xtol=1e-15, rtol=8.9e-16)
be = beta_of(s)
L  = s/lam

# hand-derived equivalent: cosh s cos s + 1 = k sinh s sin s
F = np.cosh(s)*np.cos(s) + 1 - k*np.sinh(s)*np.sin(s)
print(f"s = lam*L = {s:.10f},  F(s) = {F:.3e}")
print(f"L = {1e3*L:.6f} mm")

def u(x):
    sx = lam*x
    return A*np.cosh(sx) + be*np.sinh(sx) + C*np.cos(sx) - be*np.sin(sx)

def up(x):
    sx = lam*x
    return lam*(A*np.sinh(sx) + be*np.cosh(sx) - C*np.sin(sx) - be*np.cos(sx))

# moments (per unit width)
upp0 = lam**2*(A - C)            # u''(0)/... = lam^2 (A - C) = -k*lam^2
mtip = dbar*dc*(-upp0)           # = sqrt(2*dbar*G0)
uppL = lam**2*(A*np.cosh(s) + be*np.sinh(s) - C*np.cos(s) + be*np.sin(s))
m_ap = dbar*dc*(-uppL)           # applied end moment per unit width = moment at x=L
M    = m_ap*b

print(f"m_tip  = {mtip:.8f} N  (check vs sqrt(2 dbar G0) = {np.sqrt(2*dbar*G0):.8f})")
print(f"m_appl = {m_ap:.8f} N  (check vs sqrt(2 dbar (G0+sb*dc/2)) = {np.sqrt(2*dbar*(G0+0.5*sb*dc)):.8f})")
print(f"1) M   = {M:.6f} N m   per strip")

uh   = u(L/2)
dhaf = dc*(1-uh)
print(f"3) delta(L/2) = {uh:.9f}  ->  {1e6*dhaf:.4f} um")

# ---------------- checks ----------------
# (a) J-integral balance: m^2/(2 dbar) = G0 + int_0^dc t dd = G0 + sb*dc/2
lhs = m_ap**2/(2*dbar); rhs = G0 + 0.5*sb*dc
print(f"check J: {lhs:.10f} vs {rhs:.10f}  (rel diff {abs(lhs-rhs)/rhs:.2e})")

# (b) statics: m_appl - m_tip = int_0^L t(x) x dx
Ib = quad(lambda x: sb*u(x)*x, 0.0, L, epsabs=1e-14, epsrel=1e-13)[0]
print(f"check statics: m-mtip = {m_ap-mtip:.9f},  int t x dx = {Ib:.9f}  (diff {m_ap-mtip-Ib:.2e})")

# (c) ODE satisfied, opening monotone, 0 <= delta <= dc on [0,L]
xs = np.linspace(0.0, L, 20001)
uu = u(xs); pp = up(xs)
print(f"monotone u decreasing: {np.all(np.diff(uu) <= 1e-14)}")
print(f"u in [0,1]: min {uu.min():.3e}, max {uu.max():.6f}")
print(f"max |u'(x)| = {np.abs(pp).max():.4f} (1/m)  -> delta' <= {dc*np.abs(pp).max():.4f} mm/mm")

# (d) residual of ODE at midspan by finite difference of u''''
x0 = L/2; h = 1e-6
stencil = u(np.array([x0-2*h, x0-h, x0, x0+h, x0+2*h]))
u4 = (stencil[0] - 4*stencil[1] + 6*stencil[2] - 4*stencil[3] + stencil[4])/h**4
print(f"ODE residual at L/2: {u4 - lam**4*u(x0):.3e}")

# (e) profile numbers
print(f"delta(L) = {1e3*dc*(1-u(L)):.6f} mm (= delta_c)")
print(f"traction at front t(0) = {sb*u(0):.4f} Pa, at root t(L) = {sb*u(L):.4e} Pa")
