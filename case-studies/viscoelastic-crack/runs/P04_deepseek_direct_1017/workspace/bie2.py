import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import toeplitz, solve_toeplitz

def g(u):
    u=np.asarray(u,dtype=float)
    out=np.zeros_like(u)
    m=np.abs(u)>1e-14
    out[m]=(1/np.pi**2)*np.log(np.abs(1/np.tanh(u[m]/4)))
    return out

def T_entry(d,h,nq=16):
    a=(d-0.5)*h; b=(d+0.5)*h
    if d==0:
        # analytic: (1/pi^2) int_{-h/2}^{h/2} ln|coth(u/4)| du
        # use subtraction g = (1/pi^2)ln(4/|u|) + smooth h
        def Flog(u):
            return u*np.log(4/np.abs(u))-u
        val=(Flog(b)-Flog(a))/np.pi**2
        # smooth residual h(u)=g-(1/pi^2)ln(4/|u|)
        xg,wg=leggauss(nq)
        mid=(a+b)/2; half=(b-a)/2
        us=mid+half*xg
        hfun=g(us)-(1/np.pi**2)*np.log(4/np.abs(us))
        val+=half*np.sum(wg*hfun)
        return val
    xg,wg=leggauss(nq)
    mid=(a+b)/2; half=(b-a)/2
    us=mid+half*xg
    return half*np.sum(wg*g(us))

def b_vec(X,eps,h_zone=0.001):
    # -int_{-eps}^0 g(X-s) ds, Gauss-Legendre on zone
    nz=max(20,int(eps/h_zone))
    edges=np.linspace(-eps,0,nz+1)
    xg,wg=leggauss(8)
    b=np.zeros_like(X)
    for l,r in zip(edges[:-1],edges[1:]):
        mid=(l+r)/2; half=(r-l)/2
        ss=mid+half*xg
        for i,x in enumerate(X):
            b[i]+= -half*np.sum(wg*g(x-ss))
    return b

def solve(eps,L=15.0,N=800):
    h=L/N
    X=(np.arange(N)+0.5)*h
    dvals=np.arange(-(N-1),N)
    T=np.array([T_entry(d,h) for d in dvals])
    # T[d] with d index offset N-1; toeplitz needs first column T[d for d>=0 desc], first row
    col=T[N-1:]   # d=0..N-1
    row=T[N-1::-1] # d=0,-1,... -(N-1)
    A=toeplitz(col,row)
    b=b_vec(X,eps)
    S=solve_toeplitz((col,row),b)
    return X,S,A,b

if __name__=='__main__':
    for N in [200,400,800]:
        for eps in [0.1,0.136,0.2,0.3]:
            X,S,A,b=solve(eps,L=15,N=N)
            print('N',N,'eps',eps,'S0',S[0],'SL',S[-1],'max',S.max(),'min',S.min())
