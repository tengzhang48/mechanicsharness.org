import numpy as np
import numpy.polynomial.polynomial as P
import numpy.polynomial.chebyshev as Ch
from scipy.linalg import eig
from explore_eig import clencurt
from op_eval import build_Op_eval

def basis(K):
    base=P.Polynomial([0.0,1.0,-2.0,1.0])   # x(1-x)^2
    F=[]
    for n in range(K):
        e=np.zeros(K); e[n]=1.0
        Tn=P.Polynomial(Ch.cheb2poly(e))
        F.append(base*Tn)
    return F

def solve(K, N, nq):
    x=np.sort((np.cos(np.pi*np.arange(N+1)/N)+1)/2)   # increasing Chebyshev-Lobatto
    xe,we=np.polynomial.legendre.leggauss(nq); xe=0.5*(xe+1); we=0.5*we     # evaluation / inner-product
    xq,wq=np.polynomial.legendre.leggauss(nq+40); xq=0.5*(xq+1); wq=0.5*wq  # integration
    F=basis(K)
    Phi=np.array([f(xe) for f in F])           # K x ne
    P3c=np.array([f.deriv(3)(x) for f in F])   # K x (N+1) exact at Cheb nodes
    Ope=build_Op_eval(x,N,xe,xq,wq)            # ne x (N+1)
    Tphi=Ope@P3c.T                             # ne x K
    R=np.zeros((K,K)); M=np.zeros((K,K))
    for m in range(K):
        R[m,:]=np.sum((we*Phi[m])[:,None]*(-Tphi),axis=0)
        M[m,:]=np.sum((we*Phi[m])[:,None]*Phi.T,axis=0)
    lam,vec=eig(R,M)
    idx=np.argsort(lam.real)
    return lam[idx],vec

if __name__=="__main__":
    for N in [80,160,240]:
        for K in [10,16,22,30,40]:
            lam,_=solve(K,N,nq=300)
            print("N=%d K=%d lam="%(N,K), np.round(lam.real[:6],4))