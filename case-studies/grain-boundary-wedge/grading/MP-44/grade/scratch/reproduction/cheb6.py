import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig
from cheb5 import Rval, f_k

def build(ks):
    M=len(ks); th=(np.arange(1,M+1))*np.pi/(M+1); Zc=0.5*(1-np.cos(th))
    Rm=np.zeros((M,M)); Mm=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j,(fam,k) in enumerate(ks):
            Rm[i,j]=Rval(k,fam,Z); Mm[i,j]=f_k(k,fam,Z)
    return Rm,Mm,Zc

def resid(ks, c, lam):
    Zt=np.linspace(0.05,0.95,13)
    f=np.array([sum(c[j]*f_k(kk[1],kk[0],Z) for j,kk in enumerate(ks)) for Z in Zt])
    Rf=np.array([sum(c[j]*Rval(kk[1],kk[0],Z) for j,kk in enumerate(ks)) for Z in Zt])
    return np.max(np.abs(Rf/f-lam))

for N in [8,10,12,14,16,18,22,26,30]:
    ks=[('A',k) for k in range(N)]
    Rm,Mm,Zc=build(ks); lam,vec=eig(Rm,Mm)
    cand=sorted([(l.real,i) for i,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    vals=[c[0] for c in cand[:3]]
    r=resid(ks,vec[:,cand[0][1]].real,cand[0][0]) if cand else None
    print("A N=%2d lams:"%N, np.round(vals,4), " resid=%.1e"%r if r is not None else "")