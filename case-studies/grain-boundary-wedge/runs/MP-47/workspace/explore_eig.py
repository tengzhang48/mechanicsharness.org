"""Eigenproblem: lambda phi(x) = -int_0^1 S(x,xi) phi'''(xi) dxi,
BCs phi'(1)=0, phi''(1)=0 (phi(0)=0 automatic since S(0,xi)=0).

S(x,xi)=1/(x-xi) - 1/(x+xi) - 2 xi (x-xi)/(x+xi)^3 ; R = S - 1/(x-xi).
PV: int S p dxi = int [p(xi)-p(x)]/(x-xi) dxi + p(x)*ln|x/(1-x)| + int R p dxi.
Chebyshev-Lobatto nodes + Clenshaw-Curtis quadrature.
"""
import numpy as np
from scipy.linalg import null_space, eig

def clencurt(N):
    theta = np.pi*np.arange(N+1)/N
    x = np.cos(theta)
    w = np.zeros(N+1); ii=np.arange(1,N); v=np.ones(N-1)
    if N%2==0:
        w[0]=w[N]=1.0/(N**2-1)
        for k in range(1,N//2):
            v = v - 2*np.cos(2*k*theta[ii])/(4*k*k-1)
        v = v - np.cos(N*theta[ii])/(N*N-1)
    else:
        w[0]=w[N]=1.0/N**2
        for k in range(1,(N-1)//2+1):
            v = v - 2*np.cos(2*k*theta[ii])/(4*k*k-1)
    w[ii]=2*v/N
    return w, x

def chebD(N):
    theta=np.pi*np.arange(N+1)/N; x=np.cos(theta)
    c=np.hstack(([2.],np.ones(N-1),[2.]))*(-1)**np.arange(N+1)
    X=np.tile(x,(N+1,1)).T; dX=X-X.T
    D=np.outer(c,1./c)/(dX+np.eye(N+1))
    D=D-np.diag(D.sum(axis=1))
    return D, x

def Rreg(x, xi):
    """regular part R(x,xi)=S(x,xi)-1/(x-xi) = -1/(x+xi) - 2*xi*(x-xi)/(x+xi)^3 ; x=eval, xi=int var"""
    return -1.0/(x+xi) - 2.0*xi*(x-xi)/(x+xi)**3

def build(x, D, D3, w, N):
    """Return A_int (rows interior, cols all nodes) and D1 row helpers."""
    Op=np.zeros((N-1, N+1))
    for ii,i in enumerate(range(1,N)):
        xi=x[i]
        d=x[i]-x; d[i]=1.0
        A=np.where(np.arange(N+1)==i, 0.0, w/d)
        Op[ii,:] += A
        Op[ii,:] += w*np.array([Rreg(xi, x[j]) for j in range(N+1)])
        Op[ii,:] += w[i]*D[i,:]
        s=np.sum(np.where(np.arange(N+1)==i, 0.0, w/d))
        Op[ii,i] += -s + np.log(abs(xi/(1-xi)))
    A_int = -Op @ D3
    return A_int

def solve(N):
    w0, xc = clencurt(N)
    x=(xc+1)/2.0
    w=w0/2.0
    D0,_=chebD(N)
    Dz=2.0*D0
    D1=Dz; D2=Dz@Dz; D3=Dz@Dz@Dz
    A_int=build(x,D1,D3,w,N)
    C=np.vstack([D1[0,:], D2[0,:]])
    B=null_space(C)
    A_red=A_int@B
    M=B[1:N,:]
    lam,vec=eig(A_red, M)
    idx=np.argsort(lam.real)
    print(f"--- N={N} ---")
    # eigenvalues with smallest |real|
    j2=np.argsort(np.abs(lam.real))
    print("  smallest |Re| eigenvalues:")
    for j in j2[:6]:
        phi=B@vec[:,j]
        print(f"    lam={lam[j].real:.8f} imag={lam[j].imag:+.1e} phi(0)={phi[0]:+.1e}")
    print("  largest real (least negative):")
    for j in idx[-4:]:
        print(f"    lam={lam[j].real:.6f} imag={lam[j].imag:+.1e}")
    return x,D1,B,lam,vec,idx

if __name__=="__main__":
    # validation of operator on p=1: int S dxi = T(x,1) (T(x,0)=0)
    N=100
    w0,xc=clencurt(N); x=(xc+1)/2; w=w0/2
    D0,_=chebD(N); Dz=2*D0; D3=Dz@Dz@Dz
    _,x2=clencurt(N)
    # quick: rebuild Op directly here for check
    A_int=build(x,Dz,D3,w,N)
    def Tf(x,xi): return np.log((x+xi)/abs(x-xi)) + 6*x/(x+xi) - 2*x**2/(x+xi)**2
    p=np.ones(N+1)
    I=(-A_int)@p   # A_int = -Op@D3 ; D3 p = 0, so can't check with p=1. use Op.
    print("skip p=1 (D3 kills it)")
    for N in [60,100,140,200]:
        solve(N)