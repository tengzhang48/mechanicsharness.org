import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig

def chebT(k,x): return np.cos(k*np.arccos(np.clip(x,-1.0,1.0)))
def chebU(k,x):
    x=np.clip(x,-1.0,1.0); th=np.arccos(x); s=np.sin(th)
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(s>1e-10, np.sin((k+1)*th)/np.where(s>1e-10,s,1.0), 0.0)

xth,wth=leggauss(2500)
xphi,wphi=leggauss(240)

def Rf(k,Z):
    # R[f_k](Z) = -( SP + RP ),  SP = PV int_0^1 g/(Z-Xi) dXi (analytic),
    # RP = int_0^1 Rreg g dXi  with Rreg=(Xi^2-4 Z Xi - Z^2)/(Z+Xi)^3
    x=2.0*Z-1.0
    SP = 0.0 if k==0 else -np.pi*chebU(k-1,x)
    ph=0.5*np.pi*(xphi+1.0); w=0.5*np.pi*wphi; s=np.sin(ph); c=np.cos(ph)
    # panel 1: Xi=Z sin^2 phi
    Xi1=Z*s*s; dXi1=2.0*Z*s*c
    g1=1.0/np.sqrt(Xi1*(1.0-Xi1))*chebT(k,2.0*Xi1-1.0)
    Rg1=(Xi1**2-4.0*Z*Xi1-Z**2)/(Z+Xi1)**3
    RP1=np.sum(w*Rg1*g1*dXi1)
    # panel 2: Xi=1-(1-Z) sin^2 phi
    Xi2=1.0-(1.0-Z)*s*s; dXi2=2.0*(1.0-Z)*s*c
    g2=1.0/np.sqrt(Xi2*(1.0-Xi2))*chebT(k,2.0*Xi2-1.0)
    Rg2=(Xi2**2-4.0*Z*Xi2-Z**2)/(Z+Xi2)**3
    RP2=np.sum(w*Rg2*g2*dXi2)
    return -(SP+RP1+RP2)

def theta_int(k,Z,hfun):
    thZ=np.arcsin(np.sqrt(Z)); th=0.5*thZ*(xth+1.0); w=0.5*thZ*wth; s=np.sin(th)**2
    return np.sum(w*hfun(s)*chebT(k,np.cos(2.0*th)))
def fk(k,Z):
    F2=theta_int(k,Z,lambda s:0.5*(Z-s)**2)
    I0=theta_int(k,1.0,lambda s:np.ones_like(s)); I1=theta_int(k,1.0,lambda s:(1.0-s))
    return F2+(I0-I1)*Z+(-I0/2.0)*Z**2

def lam1(N):
    M=N; th=(np.arange(1,M+1))*np.pi/(M+1); Zc=0.5*(1-np.cos(th))
    Rm=np.zeros((M,M)); Mm=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j in range(N):
            Rm[i,j]=Rf(j,Z); Mm[i,j]=fk(j,Z)
    lam,vec=eig(Rm,Mm)
    c=sorted([l.real for l in lam if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    return c

if __name__=="__main__":
    # sanity: compare Rf against brute-force PV integration for one k,Z
    for (k,Z) in [(1,0.2),(3,0.5),(2,0.05),(2,0.001)]:
        print("Rf(k=%d,Z=%.4g)=%.6f"%(k,Z,Rf(k,Z)))
    for N in [20,30,40,60]:
        c=lam1(N); print("N=%2d lams:"%N,np.round(c[:3],5))