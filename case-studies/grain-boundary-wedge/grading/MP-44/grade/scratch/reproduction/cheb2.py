import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig

def chebT(k,x):
    return np.cos(k*np.arccos(np.clip(x,-1,1)))
def chebU(k,x):
    x=np.clip(x,-1,1); th=np.arccos(x); s=np.sin(th)
    with np.errstate(divide='ignore',invalid='ignore'):
        v=np.where(s>1e-12, np.sin((k+1)*th)/np.where(s>1e-12,s,1.0), 0.0)
    return v

def gA(k,Z):
    x=2*Z-1
    return (1-x**2)**(-0.5)*chebT(k,x)
def gB(k,Z):
    x=2*Z-1
    return (1-x**2)**0.5*chebU(k,x)

xgx,wgx=leggauss(2500)
def SP(k,fam,x):
    if fam=='A':
        if k==0: return 0.0*x
        return np.pi*chebU(k-1,x)
    else:
        return -np.pi*chebT(k+1,x)
def RP(k,fam,x):
    y=xgx
    g=gA(k,(y+1)/2) if fam=='A' else gB(k,(y+1)/2)
    ker=-1.0/(x+y+2) - 2*(y+1)*(x-y)/(x+y+2)**3
    return np.sum(wgx*ker*g)

Ng=8001; Zs=np.linspace(0,1,Ng); dz=1.0/(Ng-1)
def f_k(k,fam,Z):
    gv=gA(k,Zs) if fam=='A' else gB(k,Zs)
    I0=np.concatenate([[0],np.cumsum(0.5*(gv[1:]+gv[:-1])*dz)])
    I1=np.concatenate([[0],np.cumsum(0.5*(I0[1:]+I0[:-1])*dz)])
    I2=np.concatenate([[0],np.cumsum(0.5*(I1[1:]+I1[:-1])*dz)])
    B=I0[-1]-I1[-1]; C=-I0[-1]/2
    fv=I2+B*Zs+C*Zs**2
    idx=np.clip(np.searchsorted(Zs,Z),1,Ng-1); t=(Z-Zs[idx-1])/dz
    return fv[idx-1]*(1-t)+fv[idx]*t

def Rval(k,fam,Z):
    x=2*Z-1
    return -(SP(k,fam,x)+RP(k,fam,x))

def brute_R(gfun,Z):
    xg2,wg2=leggauss(6000); xg2=0.5*(xg2+1); wg2=0.5*wg2
    gv=gfun(xg2)
    d=np.where(np.abs(Z-xg2)>1e-12, Z-xg2, 1.0)
    sing=np.sum(wg2*(gv-gfun(Z))/d) + gfun(Z)*np.log(Z/(1-Z))
    reg=np.sum(wg2*(xg2**2-4*Z*xg2-Z**2)/(Z+xg2)**3*gv)
    return -(sing+reg)

if __name__=="__main__":
    for fam,k in [('A',1),('B',2),('A',3),('B',0)]:
        for Z in [0.2,0.5,0.7]:
            print(fam,k,round(Z,2), "Rval=%.6f"%Rval(k,fam,Z), "brute=%.6f"%brute_R(lambda W:gA(k,W) if fam=='A' else gB(k,W),Z))