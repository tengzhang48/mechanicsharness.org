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
def fth(t):
    t=np.asarray(t,float); a=np.abs(t); v=sol.sol(a)[0]; return np.where(t>=0,v,-v)
def fph(t):
    t=np.asarray(t,float); return sol.sol(np.abs(t))[1]

# I_J with angular g (A pulled out): integrand coefficient of A^(1+1/n)
th=np.linspace(-np.pi,np.pi,400001)
f,fp=fth(th),fph(th); g=np.hypot(s*f,fp)
W1=n/(n+1.0)*C*g**((n+1.0)/n)
Tr1=C*g**(1.0/n-1.0)*(s*f)
Wx1=s*f*np.cos(th)-fp*np.sin(th)
IJ=np.trapezoid(W1*np.cos(th)-Tr1*Wx1,th)
A=(J/IJ)**(n/(n+1.0))
print("I_J=%.6f  A=%.12e"%(IJ,A))

# full physical quantities
def gam(r,t):
    f,fp=fth(t),fph(t); return A*r**(s-1)*s*f, A*r**(s-1)*fp
def tau(r,t):
    gr,gt=gam(r,t); g=np.hypot(gr,gt)
    return C*g**(1.0/n-1.0)*gr, C*g**(1.0/n-1.0)*gt
def sig_xy(x,y):
    r=np.hypot(x,y); t=np.arctan2(y,x)
    tr,tt=tau(r,t); c,si=np.cos(t),np.sin(t)
    return tr*c-tt*si, tr*si+tt*c

# constitutive residual (should be ~0)
r,t=0.13,1.1; gr,gt=gam(r,t); tr,tt=tau(r,t)
tv=np.hypot(tr,tt)
pred=(gamma0/tau0)*(tv/tau0)**(n-1)*np.array([tr,tt])
print("constitutive residual rel = %.3e"%(np.hypot(*(pred-np.array([gr,gt])))/np.hypot(gr,gt)))

# Cartesian equilibrium, positive angles only
h=1e-7; res=[];mag=[]
for ang in np.linspace(0.05,np.pi-0.05,500):
    for R in (0.1,0.3,0.6):
        x,y=R*np.cos(ang),R*np.sin(ang)
        sxzp,_=sig_xy(x+h,y); sxzm,_=sig_xy(x-h,y)
        _,syyp=sig_xy(x,y+h); _,syym=sig_xy(x,y-h)
        dsxdx=(sxzp-sxzm)/(2*h); dsydy=(syyp-syym)/(2*h)
        res.append(abs(dsxdx+dsydy)); mag.append(abs(dsxdx)+abs(dsydy))
print("max equilibrium residual / max terms = %.3e"%(max(res)/max(mag)))

# face traction tau_theta(pi)
tr,tt=tau(0.2,np.pi); print("tau_theta(pi)=%.3e MPa"%(tt))

# J direct on circles
def Jnum(R):
    a=np.linspace(-np.pi,np.pi,200001)
    gr,gt=gam(R,a); tr,tt=tau(R,a)
    W=n/(n+1.0)*(tr*gr+tt*gt)
    wx=gr*np.cos(a)-gt*np.sin(a)
    return np.trapezoid((W*np.cos(a)-tr*wx)*R,a)
for R in (0.02,0.08,0.2,0.3): print("J direct r=%.2f: %.10f"%(R,Jnum(R)))

print("\nRESULTS")
for r in (0.0200,0.0800): print("  slip dw(%.4f)=%.10f mm"%(r,2*A*r**s*float(fth(np.pi))))
r0,t0=0.0200,60*np.pi/180
tr,tt=tau(r0,t0)
print("  tau_rz(0.02,60)=%.10f  tau_thz=%.10f MPa"%(tr,tt))
print("  |gamma| there=%.8f"%np.hypot(*gam(r0,t0)))
for nm,t in (("0",0.0),("90",np.pi/2),("180",np.pi)):
    print("  r(|gamma|=0.01,%sdeg)=%.10f mm"%(nm,(0.01/(A*np.hypot(s*float(fth(t)),float(fph(t)))))**(1.0/(s-1.0))))
