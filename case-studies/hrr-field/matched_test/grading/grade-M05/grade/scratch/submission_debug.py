import numpy as np
from scipy.optimize import brentq

tau0=80.0; gam0=0.002; n=3.5; a=1.0/(n+1.0)

def theta_of_u(u): return np.pi-np.arctan(u)-np.arctan(u/n)
def F_of_u(u): return (1+u*u)**(-a/2)*(u*u+n*n)**(-a*n/2)
def u_of_theta(th):
    if th<=0: return np.inf
    if th>=np.pi: return 0.0
    hi=1.0
    while theta_of_u(hi)>th: hi*=2
    return brentq(lambda u: theta_of_u(u)-th,0,hi)
def field(th):
    if th==0.0: return 0.0,0.0,a,a
    if th==np.pi:
        F=F_of_u(0.0); return 0.0,F,0.0,a*F
    u=u_of_theta(th); F=F_of_u(u); Fp=a*u*F; H=a*F*np.sqrt(1+u*u)
    return u,F,Fp,H

def J_circle(A,r):
    def integ(th):
        uu,FF,FFp,HH=field(th)
        gam=A*r**(a-1)*HH; tau=tau0*(gam/gam0)**(1.0/n)
        W=n/(n+1)*tau*gam; tr=tau*a*FF/HH
        wx=r**(a-1)*(a*FF*np.cos(th)-FFp*np.sin(th))
        return (W*np.cos(th)-tr*wx)*r
    ths=np.linspace(1e-9,np.pi-1e-9,200000)
    vals=np.array([integ(t) for t in ths])
    return 2*np.trapezoid(vals,ths)

for A in [1.0,0.004165240777881803]:
    print("A=",A,"J_circle(r=1)=",J_circle(A,1.0)," J_circle(r=0.02)=",J_circle(A,0.02))
