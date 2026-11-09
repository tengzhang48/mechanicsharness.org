import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

tau0=80.0; gamma0=0.002; n=3.5; J=0.120
p=1.0/(n+1.0); s=1.0/(n+1.0); C=tau0*gamma0**(-1.0/n)

def rhs(theta,y):
    f,fp=y; g2=(s*f)**2+fp**2; g=np.sqrt(max(g2,1e-300)); k=1.0/n-1.0
    A22=g**k+k*g**(k-2.0)*fp**2
    b=(1.0-p)*g**k*s*f+k*g**(k-2.0)*s**2*f*fp**2
    return [fp,-b/A22]

sol=solve_ivp(rhs,[0,np.pi],[0,1],method='DOP853',rtol=1e-13,atol=1e-15,dense_output=True)
# odd extension to negative theta: f(-t)=-f(t), f'(-t)=f'(t)
def fth(t):
    t=np.asarray(t,float); return np.where(t>=0, sol.sol(np.abs(t))[0], -sol.sol(np.abs(t))[0])
def fph(t):
    t=np.asarray(t,float); return sol.sol(np.abs(t))[1]
def gth(t): return np.hypot(s*fth(t),fph(t))

# I_J and A
th=np.linspace(-np.pi,np.pi,400001)
f,fp=fth(th),fph(th); g=gth(th)
W1=n/(n+1.0)*C*g**((n+1.0)/n)
Tr=C*g**(1.0/n-1.0)*(s*f)
Wx1=s*f*np.cos(th)-fp*np.sin(th)
IJ=np.trapezoid(W1*np.cos(th)-Tr*Wx1,th)
A=(J/IJ)**(n/(n+1.0))
print("I_J=%.6f  A=%.12e"%(IJ,A))

def phys(x,y):
    r=np.hypot(x,y); t=np.arctan2(y,x)
    f,fp=fth(t),fph(t); g=np.hypot(s*f,fp)
    gr=A*r**(s-1)*s*f; gt=A*r**(s-1)*fp
    tr=C*g**(1.0/n-1.0)*gr; tt=C*g**(1.0/n-1.0)*gt
    c,si=np.cos(t),np.sin(t)
    return tr*c-tt*si, tr*si+tt*c, (tr,tt)

# Cartesian equilibrium residual, positive angles only
h=1e-6; res=[];mag=[]
for ang in np.linspace(0.05,np.pi-0.05,500):
    for R in (0.1,0.3,0.6):
        x,y=R*np.cos(ang),R*np.sin(ang)
        sxzp,syzp,_=phys(x+h,y); sxzm,_,_=phys(x-h,y)
        _,syyp,_=phys(x,y+h); _,syym,_=phys(x,y-h)
        dsxdx=(sxzp-sxzm)/(2*h); dsydy=(syyp-syym)/(2*h)
        res.append(abs(dsxdx+dsydy)); mag.append(abs(dsxdx)+abs(dsydy))
print("max |equilib residual| / max |terms| = %.3e"%(max(res)/max(mag)))

# face traction
print("tau_theta(pi)=%.3e MPa (should be ~0)"%(C*gth(np.pi)**(1.0/n-1.0)*fph(np.pi)))

# constitutive residual: gamma vs tau
gr,gt=gth(0.7),None
r=0.13; t=1.1; f,fp=fth(t),fph(t); g=np.hypot(s*f,fp)
grv=A*r**(s-1)*s*f; gtv=A*r**(s-1)*fp
_,_,(tr,tt)=phys(r*np.cos(t),r*np.sin(t))
# law: gamma = (gamma0/tau0)*(tau/tau0)^(n-1) tau
tau=np.hypot(tr,tt)
pred=(gamma0/tau0)*(tau/tau0)**(n-1)*np.array([tr,tt])
print("constitutive residual rel = %.3e"%(np.hypot(*(pred-np.array([grv,gtv])))/np.hypot(grv,gtv)))

# J direct, correct A
def Jnum(R):
    a=np.linspace(-np.pi,np.pi,200001)
    f,fp=fth(a),fph(a); g=np.hypot(s*f,fp)
    gr=A*R**(s-1)*s*f; gt=A*R**(s-1)*fp
    tr=C*g**(1.0/n-1.0)*gr; tt=C*g**(1.0/n-1.0)*gt
    W=n/(n+1.0)*(tr*gr+tt*gt)
    wx=gr*np.cos(a)-gt*np.sin(a)
    return np.trapezoid((W*np.cos(a)-tr*wx)*R,a)
for R in (0.02,0.08,0.2,0.3): print("J direct r=%.2f: %.10f"%(R,Jnum(R)))

# Results
print("\nRESULTS")
for r in (0.0200,0.0800): print("  slip dw(%.4f)=%.10f mm"%(r,2*A*r**s*fth(np.pi)))
r0,t0=0.0200,60*np.pi/180
_,_,(tr,tt)=phys(r0*np.cos(t0),r0*np.sin(t0))
print("  tau_rz(0.02,60)=%.10f  tau_thz=%.10f MPa"%(tr,tt))
for nm,t in (("0",0.0),("90",np.pi/2),("180",np.pi)):
    print("  r(|gamma|=0.01,%sdeg)=%.10f mm"%(nm,(0.01/(A*gth(t)))**(1.0/(s-1.0))))
print("  f(pi)=%.10f f(pi/2)=%.10f"%(fth(np.pi),fth(np.pi/2)))
