import numpy as np
from scipy.integrate import solve_ivp

tau0=80.0; gamma0=0.002; n=3.5; J=0.120
p=1.0/(n+1.0); s=1.0/(n+1.0); C=tau0*gamma0**(-1.0/n)

def rhs(theta,y):
    f,fp=y; g2=(s*f)**2+fp**2; g=np.sqrt(max(g2,1e-300)); k=1.0/n-1.0
    A22=g**k+k*g**(k-2.0)*fp**2
    b=(1.0-p)*g**k*s*f+k*g**(k-2.0)*s**2*f*fp**2
    return [fp,-b/A22]

# integrate BOTH ways from theta=0
solu=solve_ivp(rhs,[0, np.pi],[0,1],method='DOP853',rtol=1e-13,atol=1e-15,dense_output=True)
soll=solve_ivp(rhs,[0,-np.pi],[0,1],method='DOP853',rtol=1e-13,atol=1e-15,dense_output=True)
print("f(+pi)=%.12f  f(-pi)=%.12f  sum=%.2e"%(solu.y[0,-1],soll.y[0,-1],solu.y[0,-1]+soll.y[0,-1]))
print("f'(+pi)=%.2e  f'(-pi)=%.2e"%(solu.y[1,-1],soll.y[1,-1]))

def fth(t):
    t=np.asarray(t,float)
    return np.where(t>=0, solu.sol(t)[0], soll.sol(-np.abs(t))[0])
def fph(t):
    t=np.asarray(t,float)
    return np.where(t>=0, solu.sol(t)[1], soll.sol(-np.abs(t))[1])

# I_J over full [-pi,pi]
th=np.linspace(-np.pi,np.pi,600001)
f,fp=fth(th),fph(th); g=np.hypot(s*f,fp)
W1=n/(n+1.0)*C*g**((n+1.0)/n)
Tr1=C*g**(1.0/n-1.0)*(s*f)
Wx1=s*f*np.cos(th)-fp*np.sin(th)
IJ=np.trapezoid(W1*np.cos(th)-Tr1*Wx1,th)
A=(J/IJ)**(n/(n+1.0))
print("I_J=%.6f   A=%.12e"%(IJ,A))

print("\nFINAL NUMBERS")
for r in (0.0200,0.0800):
    print("slip dw(%.4f) = %.10f mm"%(r,2*A*r**s*float(fth(np.pi))))
r0,t0=0.0200,60*np.pi/180
f,fp=fth(t0),fph(t0); g=np.hypot(s*f,fp)
# tau_r = C |gamma|^(1/n-1) gamma_r etc. with gamma_i = A r^(s-1)...
gr=A*r0**(s-1)*s*f; gt=A*r0**(s-1)*fp; G=np.hypot(gr,gt)
tr=C*G**(1.0/n-1.0)*gr; tt=C*G**(1.0/n-1.0)*gt
print("tau_rz(0.02,60deg) = %.10f MPa"%tr)
print("tau_thz(0.02,60deg) = %.10f MPa"%tt)
for nm,t in (("0",0.0),("90",np.pi/2),("180",np.pi)):
    gg=np.hypot(s*float(fth(t)),float(fph(t)))
    print("r where |gamma|=0.01 at theta=%s deg: %.10f mm"%(nm,(0.01/(A*gg))**(1.0/(s-1.0))))
