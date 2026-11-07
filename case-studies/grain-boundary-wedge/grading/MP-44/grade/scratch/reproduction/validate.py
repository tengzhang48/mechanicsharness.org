import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig
from solve3 import Rf, fk, chebT, chebU
import mpmath as mp
mp.mp.dps=25

def Rf_mp(k,Z):
    Z=mp.mpf(Z)
    N=700; x,w=leggauss(N); xi=0.5*(x+1.0); wt=0.5*w
    gv=np.array([float((mp.mpf(v)*(1-mp.mpf(v)))**(-mp.mpf(1)/2)*mp.chebyt(k,2*mp.mpf(v)-1)) for v in xi])
    gZ=float((Z*(1-Z))**(-mp.mpf(1)/2)*mp.chebyt(k,2*Z-1))
    pv=np.sum(wt*(gv-gZ)/(float(Z)-xi))+gZ*np.log(abs(float(Z)/(1-float(Z))))
    Rg=(xi**2-4*float(Z)*xi-float(Z)**2)/(float(Z)+xi)**3
    reg=np.sum(wt*Rg*gv)
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
# independent brute-force residual on f
Zt=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
print("Z     f/f(1)     R[f]/f(cheby)  R[f]/f(brute-mp)")
for Z in Zt:
    fv=sum(c[j]*fk(j,Z) for j in range(N)); Rv=sum(c[j]*Rf(j,Z) for j in range(N))
    fb=sum(c[j]*fk(j,Z) for j in range(N)); Rb=sum(c[j]*Rf_mp(j,Z) for j in range(N))
    print("%.2f  %.6f   %.6f    %.6f"%(Z,fv/f1,Rv/fv,Rb/fb))
half=sum(c[j]*fk(j,0.5) for j in range(N))/f1
print("mode1_shape_half = %.6f"%half)
print("lambda1=%.5f lambda2=%.5f"%(lam1,cand[1][0]))