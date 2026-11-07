import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig
from eig4 import build_RM

xg, wg = leggauss(80)
def Rreg(Z,X): return (X**2-4*Z*X-Z**2)/(Z+X)**3

def make_fun(N):
    edges=np.linspace(0,1,N+1); a=edges[:-1]; b=edges[1:]
    def fvals(q,Z):
        out=0.0
        for j in range(N):
            aj,bj=a[j],b[j]
            def F0(z):
                return max(0.0,min(z,bj)-aj)
            def F1(z):
                if z<=aj: return 0.0
                zz=min(z,bj); v=0.5*(zz-aj)**2
                if z>bj: v+=(bj-aj)*(z-bj)
                return v
            def F2(z):
                if z<=aj: return 0.0
                zz=min(z,bj); v=(1/6)*(zz-aj)**3
                if z>bj:
                    L=bj-aj; v+=0.5*L**2*(z-bj)+0.5*L*(z-bj)**2
                return v
            F0_1=F0(1.0);F1_1=F1(1.0);B=F0_1-F1_1;C=-F0_1/2
            out+=q[j]*(F2(Z)+B*Z+C*Z**2)
        return out
    def Rvals(q,Z):
        pv=np.log(np.abs(Z-a)/np.abs(Z-b))
        out=0.0
        for j in range(N):
            Xj=0.5*(b[j]-a[j])*xg+0.5*(a[j]+b[j]); wr=0.5*(b[j]-a[j])*wg
            reg=np.sum(wr*Rreg(Z,Xj))
            out+=q[j]*(pv[j]+reg)
        return -out
    return fvals,Rvals

for N in [160]:
    Rmat,Mmat,Zc=build_RM(N)
    lam,vec=eig(Rmat,Mmat)
    cand=sorted([(l.real,k) for k,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    fvals,Rvals=make_fun(N)
    for (lr,k) in cand[:3]:
        q=vec[:,k].real
        Zt=np.linspace(0.03,0.97,15)
        f=np.array([fvals(q,Z) for Z in Zt])
        Rf=np.array([Rvals(q,Z) for Z in Zt])
        print("lam=%.5f"%lr)
        print("  q at 15 pts:",np.round(q[::N//15][:15],4))
        print("  f:",np.round(f,5))
        print("  R[f]/f:",np.round(Rf/f,3))