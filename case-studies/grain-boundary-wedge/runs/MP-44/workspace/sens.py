import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig

def chebT(k,x): return np.cos(k*np.arccos(np.clip(x,-1.0,1.0)))
def chebU(k,x):
    x=np.clip(x,-1.0,1.0); th=np.arccos(x); s=np.sin(th)
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(s>1e-10, np.sin((k+1)*th)/np.where(s>1e-10,s,1.0), 0.0)

xth,wth=leggauss(2500)
def lambda1(N,nreg,nth=None):
    xgx,wgx=leggauss(nreg)
    def gA(k,y): return (1.0-y**2)**(-0.5)*chebT(k,y)
    def SP(k,x): return 0.0*x if k==0 else -np.pi*chebU(k-1,x)
    def RP(k,x):
        y=xgx; g=gA(k,y)
        return np.sum(wgx*(-1.0/(x+y+2.0)-2.0*(y+1.0)*(x-y)/(x+y+2.0)**3)*g)
    def Rf(k,Z):
        x=2.0*Z-1.0; return -(SP(k,x)+RP(k,x))
    def theta_int(k,Z,hfun):
        thZ=np.arcsin(np.sqrt(Z)); th=0.5*thZ*(xth+1.0); w=0.5*thZ*wth; s=np.sin(th)**2
        return np.sum(w*hfun(s)*chebT(k,np.cos(2.0*th)))
    def fk(k,Z):
        F2=theta_int(k,Z,lambda s:0.5*(Z-s)**2)
        I0=theta_int(k,1.0,lambda s:np.ones_like(s)); I1=theta_int(k,1.0,lambda s:(1.0-s))
        return F2+(I0-I1)*Z+(-I0/2.0)*Z**2
    M=N; th=(np.arange(1,M+1))*np.pi/(M+1); Zc=0.5*(1-np.cos(th))
    Rm=np.zeros((M,M)); Mm=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j in range(N):
            Rm[i,j]=Rf(j,Z); Mm[i,j]=fk(j,Z)
    lam,vec=eig(Rm,Mm)
    c=sorted([l.real for l in lam if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    return c[0] if c else np.nan

for nreg in [200,400,800,1600,3200]:
    print("nreg=%5d  lambda1(N=40)=%.6f"%(nreg,lambda1(40,nreg)))