import numpy as np
from scipy.linalg import null_space, eig
from explore_eig import clencurt, chebD, build

def lamphi(N):
    w0,xc=clencurt(N); x=(xc+1)/2.0; w=w0/2.0
    D0,_=chebD(N); Dz=2*D0; D1=Dz; D2=Dz@Dz; D3=Dz@Dz@Dz
    Aint=build(x,D1,D3,w,N)
    C=np.vstack([D1[0,:],D2[0,:]]); B=null_space(C)
    lam=eig(Aint@B, B[1:N,:], right=False)
    return x,lam

for N in [40,60,80,120,160,200,240]:
    x,lam=lamphi(N)
    r=lam.real
    sel=np.abs(r)<5000
    vals=np.sort(r[sel])
    print("N=%d |Re|<5000 ->"%N, np.round(vals,3))
    sel2=(r>0)&(r<1e6)
    print("   positive small:", np.round(np.sort(r[sel2])[:6],4))