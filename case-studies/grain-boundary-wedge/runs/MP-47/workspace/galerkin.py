"""Galerkin eigen-solver.

Basis  Phi_n(x) = x (1-x)^2 T_n(2x-1),  n=0..K-1.
These satisfy phi(0)=0, phi'(1)=0, phi''(1)=0 exactly.
Derivatives are taken on polynomials (exact, no differentiation noise).
The integral operator Op= int_0^1 S(x,xi) (...) dxi is the validated one.
Equation:  lambda phi(x) = - int S(x,xi) phi'''(xi) dxi.
Galerkin:  lambda <Phi_m,phi> = <Phi_m, -Op phi'''>   with L2 inner product.
"""
import numpy as np
from scipy.linalg import eig
import numpy.polynomial.polynomial as P
import numpy.polynomial.chebyshev as Ch
from explore_eig import clencurt, chebD, Rreg

def build_Op(x,w,D1,N):
    Op=np.zeros((N-1,N+1))
    for ii,i in enumerate(range(1,N)):
        xi=x[i]
        d=x[i]-x; d[i]=1.0
        A=np.where(np.arange(N+1)==i,0.0,w/d)
        Op[ii,:]+=A
        Op[ii,:]+=w*np.array([Rreg(xi,x[j]) for j in range(N+1)])
        Op[ii,:]+=w[i]*D1[i,:]
        s=np.sum(np.where(np.arange(N+1)==i,0.0,w/d))
        Op[ii,i]+= -s+np.log(abs(xi/(1-xi)))
    return Op

def basis(K):
    base=P.Polynomial([0.0,1.0,-2.0,1.0])   # x(1-x)^2
    F=[]
    for n in range(K):
        e=np.zeros(K); e[n]=1.0
        Tn=P.Polynomial(Ch.cheb2poly(e))
        F.append(base*Tn)
    return F

def solve(K,N):
    w0,xc=clencurt(N); x=(xc+1)/2.0; w=w0/2.0
    D0,_=chebD(N); Dz=2*D0
    Op=build_Op(x,w,Dz,N)
    F=basis(K)
    Phi=np.array([f(x) for f in F])            # K x (N+1)
    P3=np.array([f.deriv(3)(x) for f in F])    # K x (N+1)
    OPP3=Op@P3.T                                # (N-1) x K   = int S Phi_n''' at interior nodes
    R=np.zeros((K,K)); M=np.zeros((K,K))
    wi=w[1:N]                                   # interior weights
    for m in range(K):
        R[m,:]=np.sum((wi*Phi[m,1:N])[:,None]*(-OPP3),axis=0)
        M[m,:]=np.sum((wi*Phi[m,1:N])[:,None]*Phi[:,1:N].T,axis=0)
    lam,vec=eig(R,M)
    idx=np.argsort(-lam.real)
    return lam[idx],vec[:,idx]

if __name__=="__main__":
    for N in [120,200,320]:
        for K in [8,12,16,20,24,30]:
            lam,_=solve(K,N)
            print("N=%d K=%d  lam="% (N,K), np.round(lam.real[:6],4))