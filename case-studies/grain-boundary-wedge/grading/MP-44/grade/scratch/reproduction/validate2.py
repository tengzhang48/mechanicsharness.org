import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig
from solve3 import Rf, fk

Nx=800
_x,_w=leggauss(Nx)
_xi=0.5*(_x+1.0); _wt=0.5*_w
_y=2.0*_xi-1.0
def _chebT(k,y):
    T0=np.ones_like(y); T1=y.copy()
    if k==0: return T0
    if k==1: return T1
    for n in range(2,k+1):
        T0,T1=T1,2.0*y*T1-T0
    return T1
_G={}
def Rf_mp(k,Z):
    # fully independent: PV by singularity subtraction, no analytic identity
    if k not in _G: _G[k]=_chebT(k,_y)/np.sqrt(_xi*(1.0-_xi))
    gv=_G[k]
    gZ=_chebT(k,np.array([2.0*Z-1.0]))[0]/np.sqrt(Z*(1.0-Z))
    pv=np.sum(_wt*(gv-gZ)/(Z-_xi))+gZ*np.log(abs(Z/(1.0-Z)))
    Rg=(_xi**2-4.0*Z*_xi-Z**2)/(Z+_xi)**3
    reg=np.sum(_wt*Rg*gv)
    return -(pv+reg)

N=60
M=N; th=(np.arange(1,M+1))*np.pi/(M+1); Zc=0.5*(1-np.cos(th))
Rm=np.zeros((M,M)); Mm=np.zeros((M,M))
for i,Z in enumerate(Zc):
    for j in range(N):
        Rm[i,j]=Rf(j,Z); Mm[i,j]=fk(j,Z)
lam,vec=eig(Rm,Mm)
cand=sorted([(l.real,i) for i,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
print("N=%d  smallest positive eigenvalues:"%N, np.round([v[0] for v in cand[:4]],5))
lam1=cand[0][0]; c=vec[:,cand[0][1]].real
f1=sum(c[j]*fk(j,1.0) for j in range(N))
print("Z     f/f(1)     R[f]/f(cheby)  R[f]/f(brute)")
for Z in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
    fv=sum(c[j]*fk(j,Z) for j in range(N)); Rv=sum(c[j]*Rf(j,Z) for j in range(N))
    Rb=sum(c[j]*Rf_mp(j,Z) for j in range(N))
    print("%.2f  %.6f   %.6f    %.6f"%(Z,fv/f1,Rv/fv,Rb/fv))
half=sum(c[j]*fk(j,0.5) for j in range(N))/f1
print("mode1_shape_half = %.6f"%half)
print("lambda1=%.5f lambda2=%.5f lambda3=%.5f"%(lam1,cand[1][0],cand[2][0]))