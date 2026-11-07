import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig
from eig4 import build_RM

N=160
Rmat,Mmat,Zc=build_RM(N)
lam,vec=eig(Rmat,Mmat)
cand=sorted([(l.real,k) for k,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
lr,k=cand[0]; q=vec[:,k].real
edges=np.linspace(0,1,N+1)
print("lam1=",lr)
# print q near left end and right end
print("q (midpoints, first 12):", np.round(q[:12],5))
print("q (midpoints, last 12):", np.round(q[-12:],5))
print("q max abs:", np.max(np.abs(q)))
# normalize q
q=q/np.max(np.abs(q))
# check q vs Z^beta near 0
Zmid=0.5*(edges[:-1]+edges[1:])
for beta in [0.5,1.0,1.5,2.0,2.5,3.0]:
    # fit slope in log-log for small Z
    m=(Zmid>1e-3)&(Zmid<0.05)
    if m.sum()>3 and np.all(q[m]!=0):
        A=np.polyfit(np.log(Zmid[m]),np.log(np.abs(q[m])+1e-30),1)
        print("beta fit (left)=",A[0])
        break
# right end
for _ in [0]:
    m=(Zmid<0.999)&(Zmid>0.95)
    if m.sum()>3 and np.all(q[m]!=0):
        A=np.polyfit(np.log(1-Zmid[m]),np.log(np.abs(q[m])+1e-30),1)
        print("beta fit (right)=",A[0])