import numpy as np
from numpy.polynomial.legendre import leggauss

def K(u):
    u=np.asarray(u,dtype=float)
    out=np.zeros_like(u)
    m=np.abs(u)>1e-12
    out[m]=1/np.sinh(u[m]/2)
    return out

def Iconst(X,eps):
    # PV integral_{-eps}^{inf} ds/sinh((X-s)/2), finite part
    # antiderivative: -2 ln|tanh((X-s)/2)|
    # at inf -> 0 ; at -eps: -2 ln|tanh((X+eps)/2)|
    return 2*np.log(np.abs(np.tanh((X+eps)/2)))

def PVpos(X,L):
    # PV integral_0^inf ds/sinh((X-s)/2) = -2 ln|tanh(X/2)|
    return -2*np.log(np.abs(np.tanh(X/2)))

def solve(eps,L=15.0,N=600):
    X=np.linspace(0,L,N)
    w=np.ones(N); w[0]=0.5; w[-1]=0.5; w*=L/(N-1)
    A=np.zeros((N-1,N-1))
    for i in range(1,N):
        u=X[i]-X
        row=w*K(u)
        row[i]=0.0
        diag=PVpos(X[i],L)-np.sum(row)
        row[i]=diag
        A[i-1,:]=row[1:]
    b=-Iconst(X[1:],eps)
    St=np.linalg.solve(A,b)
    S=np.concatenate([[0.0],St])
    return X,S

if __name__=='__main__':
    for N in [200,400,800,1600]:
        for eps in [0.1,0.136,0.2,0.3]:
            X,S=solve(eps,L=15,N=N)
            print('N',N,'eps',eps,'Stilde0',S[0],'StildeL',S[-1],'S0',1+S[0],'SL',1+S[-1],'min',(1+S).min())
