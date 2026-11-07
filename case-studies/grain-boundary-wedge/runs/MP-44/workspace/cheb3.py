import numpy as np, time
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig
tp=time.perf_counter()

def chebT(k,x): return np.cos(k*np.arccos(np.clip(x,-1,1)))
def chebU(k,x):
    x=np.clip(x,-1,1); th=np.arccos(x); s=np.sin(th)
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(s>1e-10, np.sin((k+1)*th)/np.where(s>1e-10,s,1.0), 0.0)

def gA(k,Z):
    x=2*Z-1; return (1-x**2)**(-0.5)*chebT(k,x)
def gB(k,Z):
    x=2*Z-1; return (1-x**2)**0.5*chebU(k,x)

xgx,wgx=leggauss(600)
def SP(k,fam,x):
    if fam=='A':
        return 0.0*x if k==0 else np.pi*chebU(k-1,x)
    return -np.pi*chebT(k+1,x)
def RP(k,fam,x):
    y=xgx; g=gA(k,(y+1)/2) if fam=='A' else gB(k,(y+1)/2)
    return np.sum(wgx*(-1.0/(x+y+2) - 2*(y+1)*(x-y)/(x+y+2)**3)*g)
def Rval(k,fam,Z):
    x=2*Z-1; return -(SP(k,fam,x)+RP(k,fam,x))

Ng=4001; Zs=np.linspace(0,1,Ng); dz=1.0/(Ng-1)
def f_k(k,fam,Z):
    gv=gA(k,Zs) if fam=='A' else gB(k,Zs)
    I0=np.concatenate([[0],np.cumsum(0.5*(gv[1:]+gv[:-1])*dz)])
    I1=np.concatenate([[0],np.cumsum(0.5*(I0[1:]+I0[:-1])*dz)])
    I2=np.concatenate([[0],np.cumsum(0.5*(I1[1:]+I1[:-1])*dz)])
    B=I0[-1]-I1[-1]; C=-I0[-1]/2
    fv=I2+B*Zs+C*Zs**2
    # normalize sign/magnitude
    idx=np.clip(np.searchsorted(Zs,Z),1,Ng-1); t=(Z-Zs[idx-1])/dz
    return fv[idx-1]*(1-t)+fv[idx]*t

print("setup %.1fs"%(time.perf_counter()-tp))

def brute_R(gfun,Z):
    xg2,wg2=leggauss(1500); xg2=0.5*(xg2+1); wg2=0.5*wg2
    gv=gfun(xg2)
    d=np.where(np.abs(Z-xg2)>1e-12, Z-xg2, 1.0)
    sing=np.sum(wg2*(gv-gfun(Z))/d)+gfun(Z)*np.log(Z/(1-Z))
    reg=np.sum(wg2*(xg2**2-4*Z*xg2-Z**2)/(Z+xg2)**3*gv)
    return -(sing+reg)

for fam,k in [('A',1),('B',2)]:
    for Z in [0.2,0.5]:
        print("VERIFY",fam,k,Z,"Rval=%.6f"%Rval(k,fam,Z),"brute=%.6f"%brute_R(lambda W:gA(k,W) if fam=='A' else gB(k,W),Z))

def build(ks):
    M=len(ks); th=(np.arange(1,M+1))*np.pi/(M+1); Zc=0.5*(1-np.cos(th))
    Rm=np.zeros((M,M)); Mm=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j,(fam,k) in enumerate(ks):
            Rm[i,j]=Rval(k,fam,Z); Mm[i,j]=f_k(k,fam,Z)
    return Rm,Mm,Zc

for label,ks in [("A0-9",[('A',k) for k in range(10)]),
                 ("B0-9",[('B',k) for k in range(10)]),
                 ("mix",[('A',k) for k in range(6)]+[('B',k) for k in range(6)])]:
    Rm,Mm,Zc=build(ks)
    lam,vec=eig(Rm,Mm)
    cand=sorted([(l.real,i) for i,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    print(label,"lams:",np.round([c[0] for c in cand[:4]],5), "  t=%.0fs"%(time.perf_counter()-tp))