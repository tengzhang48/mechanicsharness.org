import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig

def chebT(k,x): return np.cos(k*np.arccos(np.clip(x,-1.0,1.0)))
def chebU(k,x):
    x=np.clip(x,-1.0,1.0); th=np.arccos(x); s=np.sin(th)
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(s>1e-10, np.sin((k+1)*th)/np.where(s>1e-10,s,1.0), 0.0)

xth,wth=leggauss(3000)          # theta grid for f_k triple integrals
xphi,wphi=leggauss(300)         # phi grid for the regular kernel part

def Rf(k,Z):
    # R[f_k](Z) = -[ PV int_0^1 g/(Z-Xi) dXi + int_0^1 Rreg g dXi ]
    # PV part (exact):  PV int_0^1 (Xi(1-Xi))^{-1/2} T_k(2Xi-1)/(Z-Xi) dXi
    #                   = -2 pi U_{k-1}(x),  x=2Z-1 ;  =0 for k=0
    # Rreg(Z,Xi) = (Xi^2 - 4 Z Xi - Z^2)/(Z+Xi)^3   (the part of S~ besides 1/(Z-Xi))
    x=2.0*Z-1.0
    SP = 0.0 if k==0 else -2.0*np.pi*chebU(k-1,x)
    ph=0.25*np.pi*(xphi+1.0); w=0.25*np.pi*wphi; s=np.sin(ph); c=np.cos(ph)
    # panel 1 : Xi = Z sin^2(phi)   -> (0,Z)
    X1=Z*s*s; dX1=2.0*Z*s*c
    g1=chebT(k,2.0*X1-1.0)/np.sqrt(X1*(1.0-X1))
    RP1=np.sum(w*((X1**2-4.0*Z*X1-Z**2)/(Z+X1)**3)*g1*dX1)
    # panel 2 : Xi = 1-(1-Z) sin^2(phi) -> (Z,1)
    X2=1.0-(1.0-Z)*s*s; dX2=2.0*(1.0-Z)*s*c
    g2=chebT(k,2.0*X2-1.0)/np.sqrt(X2*(1.0-X2))
    RP2=np.sum(w*((X2**2-4.0*Z*X2-Z**2)/(Z+X2)**3)*g2*dX2)
    return -(SP+RP1+RP2)

def theta_int(k,Z,hfun):
    # int_0^Z hfun(s) g_k(s) ds ,  s=sin^2(theta);  g_k ds = 2 T_k(cos2th) dth
    thZ=np.arcsin(np.sqrt(Z)); th=0.5*thZ*(xth+1.0); w=0.5*thZ*wth; s=np.sin(th)**2
    return 2.0*np.sum(w*hfun(s)*chebT(k,np.cos(2.0*th)))
def fk(k,Z):
    # f_k = int_0^Z (Z-s)^2/2 g ds + B Z + C Z^2 , f(0)=f'(1)=f''(1)=0
    F2=theta_int(k,Z,lambda s:0.5*(Z-s)**2)
    I0=theta_int(k,1.0,lambda s:np.ones_like(s)); I1=theta_int(k,1.0,lambda s:(1.0-s))
    return F2+(I0-I1)*Z+(-I0/2.0)*Z**2

def solve(N):
    M=N; th=(np.arange(1,M+1))*np.pi/(M+1); Zc=0.5*(1-np.cos(th))
    Rm=np.zeros((M,M)); Mm=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j in range(N):
            Rm[i,j]=Rf(j,Z); Mm[i,j]=fk(j,Z)
    lam,vec=eig(Rm,Mm)
    return lam,vec,Zc

if __name__=="__main__":
    print("Rf checks (expect brute 2.296224 / 4.833361 / 5.411481 / -5.843519 / 0.163963 / -6.179436):")
    for k,Z in [(1,0.2),(1,0.5),(1,0.7),(2,0.2),(2,0.5),(3,0.5)]:
        print("  Rf(k=%d,Z=%.1f)=%.6f"%(k,Z,Rf(k,Z)))
    print("fk checks: f0(0.5)=%.6f f1(0.5)=%.6f f2(0.5)=%.6f"%(fk(0,0.5),fk(1,0.5),fk(2,0.5)))
    for N in [12,16,20,26,32,40,50,60]:
        lam,vec,Zc=solve(N)
        c=sorted([(l.real,i) for i,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
        print("N=%2d lams:"%N,np.round([v[0] for v in c[:3]],5))