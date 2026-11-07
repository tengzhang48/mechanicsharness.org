import numpy as np, time
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig
tp=time.perf_counter()

def chebT(k,x): return np.cos(k*np.arccos(np.clip(x,-1,1)))
def chebU(k,x):
    x=np.clip(x,-1,1); th=np.arccos(x); s=np.sin(th)
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(s>1e-10, np.sin((k+1)*th)/np.where(s>1e-10,s,1.0), 0.0)

# x,y quadrature on (-1,1) for regular kernel
xgx,wgx=leggauss(400)
def gA_y(k,y): return (1-y**2)**(-0.5)*chebT(k,y)
def gB_y(k,y): return (1-y**2)**0.5*chebU(k,y)

def SP(k,fam,x):
    if fam=='A':
        return 0.0*x if k==0 else -np.pi*chebU(k-1,x)
    return np.pi*chebT(k+1,x)
def RP(k,fam,x):
    y=xgx; g=gA_y(k,y) if fam=='A' else gB_y(k,y)
    return np.sum(wgx*(-1.0/(x+y+2) - 2*(y+1)*(x-y)/(x+y+2)**3)*g)
def Rval(k,fam,Z):
    x=2*Z-1; return -(SP(k,fam,x)+RP(k,fam,x))

# theta-quadrature for f_k
Nth=2000
xth,wth=leggauss(Nth)
def theta_int(k,fam,Z,hfun):
    """int_0^Z hfun(s) g_k(s) ds, s=sin^2 th, via th in [0,thZ]."""
    thZ=np.arcsin(np.sqrt(Z))
    th=0.5*thZ*(xth+1.0); w=0.5*thZ*wth
    s=np.sin(th)**2
    if fam=='A':
        return np.sum(w*hfun(s)*chebT(k,np.cos(2*th)))
    else:
        return np.sum(w*hfun(s)*chebU(k,np.cos(2*th))*np.sin(2*th)**2)
def f_k(k,fam,Z):
    I0Z=theta_int(k,fam,Z,lambda s:np.ones_like(s))
    I1Z=theta_int(k,fam,Z,lambda s:(Z-s))
    F2Z=theta_int(k,fam,Z,lambda s:0.5*(Z-s)**2)
    I0=theta_int(k,fam,1.0,lambda s:np.ones_like(s))
    I1=theta_int(k,fam,1.0,lambda s:(1.0-s))
    B=I0-I1; C=-I0/2.0
    return F2Z+B*Z+C*Z**2

print("setup %.1fs"%(time.perf_counter()-tp))

# verification
def brute_R(gfun,Z):
    xg2,wg2=leggauss(4000); xg2=0.5*(xg2+1); wg2=0.5*wg2
    gv=gfun(xg2)
    d=np.where(np.abs(Z-xg2)>1e-14, Z-xg2, 1.0)
    sing=np.sum(wg2*(gv-gfun(Z))/d)+gfun(Z)*np.log(Z/(1-Z))
    reg=np.sum(wg2*(xg2**2-4*Z*xg2-Z**2)/(Z+xg2)**3*gv)
    return -(sing+reg)
def gfun_of(fam,k):
    if fam=='A': return lambda Z: (1-(2*Z-1)**2)**(-0.5)*chebT(k,2*Z-1)
    return lambda Z: (1-(2*Z-1)**2)**0.5*chebU(k,2*Z-1)
for fam,k in [('A',1),('A',3),('B',0),('B',2)]:
    for Z in [0.2,0.5,0.7]:
        print("VERIFY",fam,k,round(Z,2),"Rval=%.6f"%Rval(k,fam,Z),"brute=%.6f"%brute_R(gfun_of(fam,k),Z))

def build(ks):
    M=len(ks); th=(np.arange(1,M+1))*np.pi/(M+1); Zc=0.5*(1-np.cos(th))
    Rm=np.zeros((M,M)); Mm=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j,(fam,k) in enumerate(ks):
            Rm[i,j]=Rval(k,fam,Z); Mm[i,j]=f_k(k,fam,Z)
    return Rm,Mm,Zc

for label,ks in [("A0-7",[('A',k) for k in range(8)]),
                 ("B0-7",[('B',k) for k in range(8)]),
                 ("mixA4B4",[('A',k) for k in range(4)]+[('B',k) for k in range(4)])]:
    Rm,Mm,Zc=build(ks)
    lam,vec=eig(Rm,Mm)
    cand=sorted([(l.real,i) for i,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    print(label,"lams:",np.round([c[0] for c in cand[:4]],5)," t=%.0fs"%(time.perf_counter()-tp))