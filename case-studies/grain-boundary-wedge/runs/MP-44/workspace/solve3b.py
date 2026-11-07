from solve3 import Rf, fk
import numpy as np
from scipy.linalg import eig
for N in [64,70,80,90]:
    M=N; th=(np.arange(1,M+1))*np.pi/(M+1); Zc=0.5*(1-np.cos(th))
    Rm=np.zeros((M,M)); Mm=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j in range(N):
            Rm[i,j]=Rf(j,Z); Mm[i,j]=fk(j,Z)
    lam,vec=eig(Rm,Mm)
    c=sorted([(l.real,i) for i,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    print("N=%2d lams:"%N,np.round([v[0] for v in c[:3]],6))