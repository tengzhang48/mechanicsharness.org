import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

# Material
tau0 = 80.0      # MPa
gam0 = 0.002
n    = 3.5
a    = 1.0/(n+1.0)      # w ~ r^a F(theta);  a = 1/(n+1)
c    = n/(n+1.0)**2

# ---- Closed-form angular solution (shape, C=1), parameter u in [0,inf) ----
# theta(u) = pi - arctan(u) - arctan(u/n)
# F(u)     = (1+u^2)^(-a/2) (u^2+n^2)^(-a n/2)
# F'(u)    = a u F
# H = sqrt(a^2 F^2 + F'^2) = a F sqrt(1+u^2)
def theta_of_u(u):
    return np.pi - np.arctan(u) - np.arctan(u/n)

def F_of_u(u):
    return (1.0+u*u)**(-a/2.0) * (u*u+n*n)**(-a*n/2.0)

def H_of_u(u):
    return a*F_of_u(u)*np.sqrt(1.0+u*u)

# robust inversion for interior theta in (0,pi)
def u_of_theta(th):
    if th <= 0.0:
        return np.inf
    if th >= np.pi:
        return 0.0
    # bracket: theta(0+)=pi, theta(inf)=0
    hi = 1.0
    while theta_of_u(hi) > th:
        hi *= 2.0
        if hi > 1e300:
            break
    return brentq(lambda u: theta_of_u(u)-th, 0.0, hi, xtol=1e-15, rtol=1e-15)

def field(th):
    if th < 0.0:
        uu,FF,FFp,HH = field(-th)
        return uu, -FF, FFp, HH          # F odd, F' even
    if th == 0.0:
        return 0.0, 0.0, a, a            # F, Fp, H ; F=0 so return u irrelevant
    if th == np.pi:
        F = F_of_u(0.0)
        return 0.0, F, 0.0, a*F
    u = u_of_theta(th)
    F = F_of_u(u); Fp = a*u*F; H = a*F*np.sqrt(1.0+u*u)
    return u, F, Fp, H

# sanity checks
print("a=%.10f  c=%.10f"%(a,c))
u,F,Fp,H = field(np.pi/2); print("theta=90: u=%.9f F=%.9f Fp=%.9f H=%.9f"%(u,F,Fp,H))
print("F(pi)=%.12f  n^{-an}=%.12f  H(0)=%.12f"%(F_of_u(0.0), n**(-a*n), a))

# ---- J via u-parameterization (A=1) ----
def integrand_u(u):
    F = F_of_u(u); Fp = a*u*F; H = a*F*np.sqrt(1.0+u*u)
    r = 1.0
    gam = r**(a-1.0)*H
    tau = tau0*(gam/gam0)**(1.0/n)
    W = n/(n+1.0)*tau*gam
    tr = tau*a*F/H
    den = np.sqrt((1.0+u*u)*(1.0+u*u/(n*n)))
    sinpsi = u*(1.0+1.0/n)/den
    cospsi = (1.0-u*u/n)/den
    sinth = sinpsi          # sin(pi-psi)=sin psi
    costh = -cospsi         # cos(pi-psi)=-cos psi
    wx = r**(a-1.0)*(a*F*costh - Fp*sinth)
    # -dtheta/du
    dtheta = 1.0/(1.0+u*u) + n/(n*n+u*u)
    return (W*costh - tr*wx)*r * dtheta   # times |dtheta/du|, A=1

# J = 2 * integral_0^pi integrand dtheta = 2 * integral_0^inf integrand_u du
Dfull,_ = quad(integrand_u, 0.0, np.inf, limit=600, epsabs=1e-13, epsrel=1e-13)
J_per_A = 2.0*Dfull     # includes tau0
print("J (A=1) =", J_per_A)

Jval = 0.120
A = (Jval/J_per_A)**(n/(n+1.0))
print("A =", A)

# ---- Requested quantities ----
Fpi = F_of_u(0.0)
print("\n1) face slip")
for r in [0.0200, 0.0800]:
    dw = 2.0*A*r**a*Fpi
    print("   Delta w(r=%.4f) = %.10g mm"%(r,dw))

print("\n2) stresses at r=0.02, theta=60 deg")
r=0.0200; th=np.deg2rad(60.0)
u,F,Fp,H=field(th)
gam=A*r**(a-1.0)*H; tau=tau0*(gam/gam0)**(1.0/n)
trz=tau*a*F/H; ttz=tau*Fp/H
print("   tau_rz  = %.10g MPa"%trz)
print("   tau_thz = %.10g MPa"%ttz)

print("\n3) radii where |gamma|=0.0100")
g_target=0.0100
for name,th in [("theta=0",0.0),("theta=90",np.pi/2),("theta=180",np.pi)]:
    uu,F,Fp,H=field(th)
    r_g=(A*H/g_target)**((n+1.0)/n)
    print("   %s : r = %.10g mm  (H=%.10g)"%(name,r_g,H))

# ---- Verifications ----
print("\nVerification: path independence of J on circles")
def J_circle(r):
    def integ(th):
        uu,FF,FFp,HH=field(th)
        gam=A*r**(a-1.0)*HH; tau=tau0*(gam/gam0)**(1.0/n)
        W=n/(n+1.0)*tau*gam; tr=tau*a*FF/HH
        wx=A*r**(a-1.0)*(a*FF*np.cos(th)-FFp*np.sin(th))
        return (W*np.cos(th)-tr*wx)*r
    return 2.0*quad(integ,0.0,np.pi,limit=400,epsabs=1e-12,epsrel=1e-12)[0]

for r in [0.005,0.02,0.1,0.3]:
    print("   r=%.3f  J=%.9g"%(r,J_circle(r)))

print("\nVerification: equilibrium residual div(tau) (should be ~0)")
def tau_vec(x,y):
    r=np.hypot(x,y); th=np.arctan2(y,x)
    uu,FF,FFp,HH=field(th)
    gam=A*r**(a-1.0)*HH; tau=tau0*(gam/gam0)**(1.0/n)
    er=np.array([np.cos(th),np.sin(th)]); eth=np.array([-np.sin(th),np.cos(th)])
    return tau/HH*(a*FF*er+FFp*eth)
h=1e-6
for (x,y) in [(0.02,0.02),(0.02,-0.03),(-0.05,0.04),(0.1,0.1),(0.01,0.0)]:
    r=np.hypot(x,y)
    dvx=(tau_vec(x+h,y)[0]-tau_vec(x-h,y)[0])/(2*h)
    dvy=(tau_vec(x,y+h)[1]-tau_vec(x,y-h)[1])/(2*h)
    print("   (%.2f,%.2f) div=%.3e  tau_mag=%.4g"%(x,y,dvx+dvy,np.hypot(*tau_vec(x,y))))

print("\nVerification: traction tau_thz on face theta=pi (limit)")
# F'->0 at face
print("   Fp at theta=pi =", field(np.pi)[2])

print("\nVerification: integrate angular ODE directly (scipy) and compare to closed form")
from scipy.integrate import solve_ivp
p_exp=(n+1.0)/n
def ode(th,Y):
    F,Fp=Y
    num = -F*(c*a*a*F*F + (c+(p_exp-2.0)*a*a)*Fp*Fp)
    den = a*a*F*F + (p_exp-1.0)*Fp*Fp
    return [Fp, num/den]
sol=solve_ivp(ode,[0,np.pi],[0.0,1.0],rtol=1e-12,atol=1e-14,dense_output=True,max_step=1e-3)
print("   IVP F(pi)=%.10f F'(pi)=%.3e"%(sol.y[0,-1],sol.y[1,-1]))
print("   closed form F(pi)=%.10f  F(0+)=0 F'(0+)=%.10f"%(F_of_u(0.0),a))
for th in [np.pi/3,np.pi/2,2*np.pi/3]:
    _,F,Fp,H=field(th)
    Fi=sol.sol(th)[0]; Fpi_=sol.sol(th)[1]
    print("   th=%.1f deg  F_closed=%.9f F_ivp=%.9f  Fp_closed=%.9f Fp_ivp=%.9f"%(np.rad2deg(th),F,Fi,Fp,Fpi_))

print("\nVerification: stored-work law W = n/(n+1) tau*gamma  (integrate dW=tau dgamma)")
from scipy.integrate import quad as _q
g=0.004
_tau=lambda gg: tau0*(gg/gam0)**(1.0/n)
Wnum=_q(lambda gg: _tau(gg),0,g)[0]
Wform=n/(n+1.0)*_tau(g)*g
print("   W numeric=%.12g  W formula=%.12g"%(Wnum,Wform))
