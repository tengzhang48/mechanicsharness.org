import numpy as np, time
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig

def chebT(k,x): return np.cos(k*np.arccos(np.clip(x,-1,1)))
def chebU(k,x):
    x=np.clip(x,-1,1); th=np.arccos(x); s=np.sin(th)
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(s>1e-10, np.sin((k+1)*th)/np.where(s>1e-10,s,1.0), 0.0)
xgx,wgx=leggauss(300)
def gA_y(k,y): return (1-y**2)**(-0.5)*chebT(k,y)
def gB_y(k,y): return (1-y**2)**0.5*chebU(k,y)
def SP(k,fam,x):
    if fam=='A': return 0.0*x if k==0 else -np.pi*chebU(k-1,x)
    return np.pi*chebT(k+1,x)
def RP(k,fam,x):
    y=xgx; g=gA_y(k,y) if fam=='A' else gB_y(k,y)
    return np.sum(wgx*(-1.0/(x+y+2)-2*(y+1)*(x-y)/(x+y+2)**3)*g)
def Rval(k,fam,Z):
    x=2*Z-1; return -(SP(k,fam,x)+RP(k,fam,x))
xth,wth=leggauss(2000)
def theta_int(k,fam,Z,hfun):
    thZ=np.arcsin(np.sqrt(Z)); th=0.5*thZ*(xth+1.0); w=0.5*thZ*wth; s=np.sin(th)**2
    if fam=='A': return np.sum(w*hfun(s)*chebT(k,np.cos(2*th)))
    return np.sum(w*hfun(s)*chebU(k,np.cos(2*th))*np.sin(2*th)**2)
def f_k(k,fam,Z):
    F2=theta_int(k,fam,Z,lambda s:0.5*(Z-s)**2)
    I0=theta_int(k,fam,1.0,lambda s:np.ones_like(s)); I1=theta_int(k,fam,1.0,lambda s:(1.0-s))
    return F2+(I0-I1)*Z+(-I0/2.0)*Z**2
def build(ks):
    M=len(ks); th=(np.arange(1,M+1))*np.pi/(M+1); Zc=0.5*(1-np.cos(th))
    Rm=np.zeros((M,M)); Mm=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j,(fam,k) in enumerate(ks):
            Rm[i,j]=Rval(k,fam,Z); Mm[i,j]=f_k(k,fam,Z)
    return Rm,Mm,Zc
def check(label,ks):
    Rm,Mm,Zc=build(ks); lam,vec=eig(Rm,Mm)
    cand=sorted([(l.real,i) for i,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    print(label,"lams:",np.round([c[0] for c in cand[:4]],5))
    if cand:
        lr,i=cand[0]; c=vec[:,i].real
        Zt=np.linspace(0.05,0.95,11)
        f=np.array([sum(c[j]*f_k(kk[1],kk[0],Z) for j,kk in enumerate(ks)) for Z in Zt])
        Rf=np.array([sum(c[j]*Rval(kk[1],kk[0],Z) for j,kk in enumerate(ks)) for Z in Zt])
        f1=sum(c[j]*f_k(kk[1],kk[0],1.0) for j,kk in enumerate(ks))
        print("  f(0.05..0.95):",np.round(f/f1,4))
        print("  R[f]/f:",np.round(Rf/f,3))
if __name__=="__main__":
    for N in [8,10,12]:
        check("A0-%d"%(N-1),[('A',k) for k in range(N)])
    for N in [8,10,12]:
        check("B0-%d"%(N-1),[('B',k) for k in range(N)])