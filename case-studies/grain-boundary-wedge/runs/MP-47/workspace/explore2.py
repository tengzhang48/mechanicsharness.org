import numpy as np
from scipy.linalg import null_space, eig
from explore_eig import clencurt, chebD, Rreg, build

def opmat(N):
    w0,xc=clencurt(N); x=(xc+1)/2.0; w=w0/2.0
    D0,_=chebD(N); Dz=2*D0
    D1=Dz; D2=Dz@Dz; D3=Dz@Dz@Dz
    Op=build(x,D1,D3,w,N)  # = -Op@D3 (rows interior)
    return x,D1,D2,D3,Op

def eig_phi(N):
    x,D1,D2,D3,Aint=opmat(N)
    C=np.vstack([D1[0,:],D2[0,:]])
    B=null_space(C)
    Ared=Aint@B
    M=B[1:N,:]
    lam=eig(Ared,M,right=False)
    return x,lam

def eig_q(N):
    x,D1,D2,D3,Aint=opmat(N)
    Op=-Aint@np.linalg.inv(D3)   # recover Op = int S (rows interior), since Aint=-Op@D3
    D3inv=np.linalg.inv(D3)
    C=np.vstack([D1[0,:]@D3inv, D2[0,:]@D3inv])
    Pc=null_space(C)
    A=Op@Pc
    Phil=D3inv@Pc
    lam=eig(A,Phil[1:N,:],right=False)
    return x,lam

for N in [40,60,80,120,160,200]:
    x,lp=eig_phi(N)
    x,lq=eig_q(N)
    sp=np.sort(lp.real); sq=np.sort(lq.real)
    negp=sp[sp<0]; posp=sp[sp>0]
    print(f"N={N}")
    print("  phi-form  most-neg real:", np.round(negp[-8:],4))
    print("  q-form    most-neg real:", np.round(np.sort(lq.real)[np.sort(lq.real)<0][-8:],4))