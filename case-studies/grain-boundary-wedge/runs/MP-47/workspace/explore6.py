import numpy as np
from scipy.linalg import null_space, eig
from explore_eig import clencurt, chebD
from op2 import build_Op2

def run(N, drop, ngl=600):
    w0,xc=clencurt(N); x=(xc+1)/2.0
    D0,_=chebD(N); Dz=2*D0; D1=Dz; D2=Dz@Dz; D3=Dz@Dz@Dz
    Op=build_Op2(x,N,ngl)
    Aint=-Op@D3
    eN=np.zeros(N+1); eN[N]=1.0
    C=np.vstack([eN, D1[0,:], D2[0,:]])
    B=null_space(C)
    A=Aint@B
    Mall=B[1:N,:]
    if drop=='start':   keep=np.arange(1,N-1)
    elif drop=='mid':
        keep=np.concatenate([np.arange(0,(N-2)//2),np.arange((N-2)//2+1,N-1)])
    lam=eig(A[keep,:], Mall[keep,:], right=False)
    return x,lam

for drop in ['start','mid']:
    print("drop=",drop,flush=True)
    for N in [500,650,800]:
        x,lam=run(N,drop,ngl=700)
        pos=np.sort(lam.real[(lam.real>0)])
        print("  N=%d  smallest positive:"%N, np.round(pos[:6],4),flush=True)