import numpy as np
from scipy.linalg import null_space
from explore_eig import clencurt, chebD
from op2 import build_Op2

def pencil(N, ngl=500):
    w0,xc=clencurt(N); x=(xc+1)/2.0
    D0,_=chebD(N); Dz=2*D0; D1=Dz; D2=Dz@Dz; D3=Dz@Dz@Dz
    Op=build_Op2(x,N,ngl)
    Aint=-Op@D3
    eN=np.zeros(N+1); eN[N]=1.0
    C=np.vstack([eN, D1[0,:], D2[0,:]])
    B=null_space(C)
    A=Aint@B
    M=B[1:N,:]
    return x,A,M

def smin(lam,A,M):
    s=np.linalg.svd(A-lam*M,compute_uv=False)
    return s[-1]

def scan(N, lo, hi, n):
    x,A,M=pencil(N)
    grid=np.linspace(lo,hi,n)
    vals=np.array([smin(l,A,M) for l in grid])
    return grid,vals,A,M

def refine(A,M,a,b):
    # ternary search for minimum of smin on [a,b]
    for _ in range(60):
        m1=a+(b-a)/3; m2=b-(b-a)/3
        if smin(m1,A,M)<smin(m2,A,M): b=m2
        else: a=m1
    return 0.5*(a+b), smin(0.5*(a+b),A,M)

if __name__=="__main__":
    import sys
    N=int(sys.argv[1]) if len(sys.argv)>1 else 160
    grid,vals,A,M=scan(N,0.0,420.0,20000)
    # local minima
    roots=[]
    for i in range(2,len(grid)-2):
        if vals[i]<=vals[i-1] and vals[i]<=vals[i+1] and vals[i]<0.05:
            r,s=refine(A,M,grid[i-1],grid[i+1])
            roots.append((r,s))
    # dedupe
    out=[]
    for r,s in roots:
        if not out or abs(r-out[-1][0])>1e-2:
            out.append((r,s))
    print("N=%d"%N)
    for r,s in out[:10]:
        print("   lam=%.6f  smin=%.2e"%(r,s))