import numpy as np
from scan2 import pencil

if __name__=="__main__":
    import sys
    N=int(sys.argv[1]) if len(sys.argv)>1 else 200
    x,A,M=pencil(N,ngl=600)
    grid=np.linspace(0.0,420.0,8401)          # step 0.05
    vals=np.zeros_like(grid)
    for k,l in enumerate(grid):
        vals[k]=np.linalg.svd(A-l*M,compute_uv=False)[-1]
    mins=[]
    for i in range(2,len(grid)-2):
        if vals[i]<=vals[i-1] and vals[i]<=vals[i+1]:
            mins.append((vals[i],grid[i]))
    mins.sort()
    print("N=%d  smallest local minima of sigma_min(A-lam M):"%N)
    for v,l in mins[:12]:
        print("   smin=%.3e  at lam=%.4f"%(v,l))