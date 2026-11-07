import numpy as np
from numpy.polynomial.chebyshev import chebval, chebder
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig

# work on Z in (0,1); x=2Z-1 in (-1,1)
Ng=4001
Zs=np.linspace(0,1,Ng); dz=1.0/(Ng-1)

def chebT(k,x): return np.cos(k*np.arccos(np.clip(x,-1,1)))
def chebU(k,x):
    x=np.clip(x,-1,1); th=np.arccos(x)
    return np.sin((k+1)*th)/np.sin(th)

# basis family A: g_k(Z)= 2*(1-x^2)^{-1/2} T_k(x)  (x=2Z-1)
#   (Z(1-Z))^{-1/2} = 2 (1-x^2)^{-1/2}
def gA(k,Z):
    x=2*Z-1
    w=2.0/np.sqrt(np.clip(1-x**2,1e-300,None))
    return w*chebT(k,x)
# family B: g_k(Z)= (1/2)*sqrt(1-x^2) U_k(x)   since sqrt(Z(1-Z))=sqrt(1-x^2)/2
def gB(k,Z):
    x=2*Z-1
    w=0.5*np.sqrt(np.clip(1-x**2,0,None))
    return w*chebU(k,x)

xg,wg=leggauss(3000); xg=0.5*(xg+1); wg=0.5*wg
xgp=2*xg-1

def Rreg(Z,X):
    return (X**2-4*Z*X-Z**2)/(Z+X)**3

def pv_A(k,x):
    # PV int_0^1 gA_k(Xi)/(Z-Xi) dXi ;  x=2Z-1
    if k==0: return 0.0*x
    return -np.pi*chebU(k-1,x)
def pv_B(k,x):
    return np.pi*chebT(k+1,x)

def Rfun(k,fam,Z):
    x=2*Z-1
    if fam=='A':
        pv=pv_A(k,x); g=lambda X: gA(k,X)
    else:
        pv=pv_B(k,x); g=lambda X: gB(k,X)
    reg=np.sum(wg*Rreg(Z,xg)*g(xg))
    return -(pv+reg)

def f_fun(k,fam,Z):
    gv = gA(k,Zs) if fam=='A' else gB(k,Zs)
    # cumulative integrals
    I0=np.concatenate([[0],np.cumsum(0.5*(gv[1:]+gv[:-1])*dz)])
    I1=np.concatenate([[0],np.cumsum(0.5*(I0[1:]+I0[:-1])*dz)])
    I2=np.concatenate([[0],np.cumsum(0.5*(I1[1:]+I1[:-1])*dz)])
    B=I0[-1]-I1[-1]; C=-I0[-1]/2
    fv=I2+B*Zs+C*Zs**2
    # evaluate at Z by interpolation
    idx=np.clip(np.searchsorted(Zs,Z),1,Ng-1)
    t=(Z-Zs[idx-1])/dz
    return fv[idx-1]*(1-t)+fv[idx]*t

# build pencil with basis union A0..A_NA, B0..B_NB
def build(NA,NB):
    ks=[('A',k) for k in range(NA)]+[('B',k) for k in range(NB)]
    M=len(ks)
    # collocation points (Chebyshev interior)
    th=(np.arange(1,M+1))*np.pi/(M+1); Zc=0.5*(1-np.cos(th))
    Rm=np.zeros((M,M)); Mm=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j,(fam,k) in enumerate(ks):
            Rm[i,j]=Rfun(k,fam,Z); Mm[i,j]=f_fun(k,fam,Z)
    return Rm,Mm,Zc,ks

for (NA,NB) in [(10,0),(0,10),(8,8),(12,4),(6,12),(14,14)]:
    Rm,Mm,Zc,ks=build(NA,NB)
    lam,vec=eig(Rm,Mm)
    cand=sorted([(l.real,k) for k,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    print("NA=%d NB=%d lams:"%(NA,NB), np.round([c[0] for c in cand[:4]],5))