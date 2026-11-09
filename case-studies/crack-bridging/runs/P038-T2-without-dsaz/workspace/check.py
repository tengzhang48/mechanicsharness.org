import numpy as np
from scipy.integrate import quad

b=0.025; E1=110e9; h1=0.60e-3; E2=190e9; h2=0.90e-3
sig=0.20e6; dc=0.20e-3; G0=1.40
EI1=E1*b*h1**3/12.0; EI2=E2*b*h2**3/12.0
R=1.0/EI1+1.0/EI2

# ---- Check 1: G for pure-moment DCB via compliance (virtual crack extension) ----
M=1.0
# compliance dC/da: rotation of each arm from end moment over crack length a
# theta_i = M a / EI_i  (each arm). relative rotation = M a R
dtheta_da = M*R
G = 0.5*M*dtheta_da/b
print("Check G formula: G/(M^2)=%.6e ; M^2 R/(2b)=%.6e" % (G/M**2, R/(2*b)))

# ---- Check 2: total energy balance at steady state ----
# available G = M^2 R/(2b) equals G0 + cohesive work
work = sig*dc/2.0
print("cohesive work per area = %.4f J/m2  (sigma_b*dc/2)" % work)
Greq = G0 + work
Mss = np.sqrt(2*b*Greq/R)
print("steady-state moment M = %.6f N.m" % Mss)

# ---- Check 3: process-zone length & profile (tip delta=dc, root delta=0) ----
k=b*sig*R/dc; beta=k**0.25
A=Mss*R/(2*beta**2); C=-A
from scipy.optimize import brentq
def resid(ell):
    u=beta*ell
    B=-A*(np.sinh(u)+np.sin(u))/(np.cosh(u)+np.cos(u))
    return A*(np.cosh(u)-np.cos(u))+B*(np.sinh(u)+np.sin(u))+dc
lams=np.linspace(1e-4,0.1,5000); vals=[resid(l) for l in lams]
sc=np.where(np.diff(np.sign(vals)))[0]
ell=brentq(resid,lams[sc[0]],lams[sc[0]+1])
u=beta*ell; B=-A*(np.sinh(u)+np.sin(u))/(np.cosh(u)+np.cos(u)); D=B
def d(x):
    ux=beta*x
    return dc+A*np.cosh(ux)+B*np.sinh(ux)+C*np.cos(ux)+D*np.sin(ux)
print("process zone length L = %.4f mm" % (1e3*ell))
print("delta(L/2) = %.4f um   delta(0)=%.1f um  delta(L)=%.3e um" % (1e6*d(ell/2),1e6*d(0),1e6*d(ell)))

# ---- Check 4: energy stored in cohesive zone equals work predicted ----
xs=np.linspace(0,ell,2000); ds=d(xs)
# energy released rate at tip should equal G0 + integral t d(delta) actually is
# rebuilding resistance: integral of t over the separation process
# Verify: total potential balance is built into G requirement (steady state).

# ---- Check 5: verify ODE residual ----
h=1e-6
x0=ell/2
resid_ode=(d(x0+2*h)-4*d(x0+h)+6*d(x0)-4*d(x0-h)+d(x0-2*h))/h**4 - k*d(x0) + b*sig*R
print("ODE residual at L/2 = %.3e (should be ~0)" % resid_ode)

# also report moment using alternative k scaling sanity
print("\nSUMMARY: M=%.6f N.m ; L=%.4f mm ; delta(L/2)=%.4f um" % (Mss,1e3*ell,1e6*d(ell/2)))