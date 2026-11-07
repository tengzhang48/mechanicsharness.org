import numpy as np
from scipy.linalg import null_space, eig
from explore_eig import clencurt, chebD, build

def run(N, drop):
    w0,xc=clencurt(N); x=(xc+1)/2.0; w=w0/2.0
    D0,_=chebD(N); Dz=2*D0; D1=Dz; D2=Dz@Dz; D3=Dz@Dz@Dz
    Aint=build(x,D1,D3,w,N)      # (N-1) x (N+1)
    eN=np.zeros(N+1); eN[N]=1.0
    C=np.vstack([eN, D1[0,:], D2[0,:]])
    B=null_space(C)              # (N+1) x (N-2)
    A=Aint@B                     # (N-1) x (N-2)
    Mall=B[1:N,:]                # nodes 1..N-1
    # drop 'drop' interior rows -> choose rows r0..r0+N-3
    rows=np.arange(1,N-1)        # nodes 1..N-1 (positions 0..N-2 in A)
    keep=np.arange(0,N-2)        # keep first N-2 interior rows (nodes 1..N-2)
    if drop=='end':   keep=np.arange(0,N-2)
    if drop=='start': keep=np.arange(1,N-1)
    lam=eig(A[keep,:], Mall[keep,:], right=False)
    return x,lam

for drop in ['end','start']:
    print("drop=",drop)
    for N in [50,80,120,160,200,260,320]:
        x,lam=run(N,drop)
        r=lam.real
        pos=np.sort(r[(r>0)])
        print("  N=%d  smallest positive:"%N, np.round(pos[:6],4))