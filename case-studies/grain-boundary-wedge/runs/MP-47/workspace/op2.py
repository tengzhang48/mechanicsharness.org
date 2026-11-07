import numpy as np
from explore_eig import clencurt, chebD, Rreg

def cheb_interp_matrix(x, xq, N):
    """Matrix mapping values at Chebyshev-Lobatto nodes x (N+1) to values at xq."""
    theta=np.pi*np.arange(N+1)/N
    j=np.arange(N+1)[None,:]; k=np.arange(N+1)[:,None]
    Acoef=(2.0/N)*np.cos(k*theta)      # [k,j] = (2/N) cos(k theta_j)
    Acoef[:,0]*=0.5; Acoef[:,N]*=0.5   # half weight at endpoint nodes
    t=2*xq-1.0
    T=np.cos(np.arccos(np.clip(t,-1,1))[:,None]*np.arange(N+1)[None,:])
    T[:,0]*=0.5; T[:,N]*=0.5          # double-prime: halve k=0 and k=N terms
    return T@Acoef  # (len(xq), N+1)

def build_Op2(x, N, ngl=400):
    xq,wq=np.polynomial.legendre.leggauss(ngl); xq=0.5*(xq+1); wq=0.5*wq
    Icall=cheb_interp_matrix(x,xq,N)
    Op=np.zeros((N-1,N+1))
    for ii,i in enumerate(range(1,N)):
        xi=x[i]
        coef1=wq/(xi-xq)
        Op[ii,:]+= coef1 @ Icall
        s=np.sum(coef1)
        coef3=wq*Rreg(xi,xq)
        Op[ii,:]+= coef3 @ Icall
        Op[ii,i]+= -s + np.log(abs(xi/(1.0-xi)))
    return Op

if __name__=="__main__":
    from test_op import direct
    N=80
    w0,xc=clencurt(N); x=(xc+1)/2.0
    Op=build_Op2(x,N,ngl=600)
    for pfun in [lambda t:t**2, lambda t:np.exp(2*t), lambda t:np.sin(3*t)]:
        I=Op@pfun(x); Idir=np.array([direct(xi,pfun) for xi in x[1:-1]])
        print("rel err", np.max(np.abs(I-Idir))/np.max(np.abs(Idir)))