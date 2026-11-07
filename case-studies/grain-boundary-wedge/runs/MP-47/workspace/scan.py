import numpy as np
from scipy.linalg import null_space
from explore_eig import clencurt, chebD, build

def pencil(N):
    w0,xc=clencurt(N); x=(xc+1)/2.0; w=w0/2.0
    D0,_=chebD(N); Dz=2*D0; D1=Dz; D2=Dz@Dz; D3=Dz@Dz@Dz
    Aint=build(x,D1,D3,w,N)          # (N-1) x (N+1)
    eN=np.zeros(N+1); eN[N]=1.0
    C=np.vstack([eN, D1[0,:], D2[0,:]])
    B=null_space(C)                  # (N+1) x (N-2)
    A=Aint@B                         # (N-1) x (N-2)
    M=B[1:N,:]                       # (N-1) x (N-2)
    return x,A,M

def smin(lam,A,M):
    T=A-lam*M
    s=np.linalg.svd(T,compute_uv=False)
    return s[-1]/np.linalg.svd(M,compute_uv=False)[0]

def find_roots(N, lo=0.0, hi=400.0, n=40000):
    x,A,M=pencil(N)
    lams=np.linspace(lo,hi,n)
    vals=np.array([smin(l,A,M) for l in lams])
    # local minima near zero
    roots=[]
    for i in range(1,n-1):
        if vals[i]<vals[i-1] and vals[i]<vals[i+1] and vals[i]<1e-3:
            # refine with golden/bisection on derivative -> just parabola
            a,b,c=lams[i-1],lams[i],lams[i+1]
            va,vb,vc=vals[i-1],vals[i],vals[i+1]
            den=(va-2*vb+vc)
            if abs(den)>0:
                dr=a-0.5*(b-a)*(vc-va)/den if den!=0 else b
            else: dr=b
            # use vertex of parabola through (a,va),(b,vb),(c,vc)
            # general quadratic fit
            dr=b+0.5*(b-a)*(va-vc)/(va-2*vb+vc) if (va-2*vb+vc)!=0 else b
            roots.append(dr)
    return roots, x, A, M

if __name__=="__main__":
    for N in [120,200,300,400]:
        r,x,A,M=find_roots(N,0,600)
        print("N=%d roots(<600):"%N, np.round(r[:8],4))