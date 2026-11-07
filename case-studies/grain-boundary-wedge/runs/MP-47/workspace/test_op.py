"""Validate the integral operator Op: I(x)=int_0^1 S(x,xi) p(xi) dxi against
direct PV quadrature (scipy.integrate.quad, cauchy weight is 1/(t-c))."""
import numpy as np
from scipy.integrate import quad
from explore_eig import clencurt, chebD, Rreg

def Sfun(x, xi):
    return 1.0/(x-xi) - 1.0/(x+xi) - 2.0*xi*(x-xi)/(x+xi)**3

def direct(xi_pt, pfun):
    # int S(xi_pt,xi) p(xi) dxi = PV int p/(xi_pt-xi) + int R p
    cauchy,_ = quad(lambda xi: pfun(xi), 0, 1, weight='cauchy', wvar=xi_pt)  # = PV int p/(xi-xi_pt)
    pv = -cauchy                                                             # PV int p/(xi_pt-xi)
    reg,_ = quad(lambda xi: Rreg(xi_pt,xi)*pfun(xi), 0, 1)
    return pv + reg

def build_Op(x,w,D1,N):
    Op=np.zeros((N-1,N+1))
    for ii,i in enumerate(range(1,N)):
        xi=x[i]; d=x[i]-x; d[i]=1.0
        A=np.where(np.arange(N+1)==i,0.0,w/d)
        Op[ii,:]+=A
        Op[ii,:]+=w*np.array([Rreg(xi,x[j]) for j in range(N+1)])
        Op[ii,:]+=w[i]*D1[i,:]
        s=np.sum(np.where(np.arange(N+1)==i,0.0,w/d))
        Op[ii,i]+= -s+np.log(abs(xi/(1-xi)))
    return Op

def main():
    N=80
    w0,xc=clencurt(N); x=(xc+1)/2; w=w0/2
    D0,_=chebD(N); Dz=2*D0
    Op=build_Op(x,w,Dz,N)
    xs=x[1:-1]
    for name,pfun,p3fun in [("p=x^2", lambda t:t**2, lambda t:0*t),
                            ("p=exp(2x)", lambda t:np.exp(2*t), lambda t:8*np.exp(2*t))]:
        pdisc=pfun(x)
        I=np.zeros(N-1)
        if name=="p=x^2":
            I=Op@pdisc      # Op acts on p directly (integral of S*p)
            Idir=np.array([direct(xi,pfun) for xi in xs])
        else:
            # integral of S * p'''(xi) : act Op on p'''
            p3d=(Dz@Dz@Dz)@pdisc
            I=Op@p3d
            Idir=np.array([direct(xi,p3fun) for xi in xs])
        err=np.max(np.abs(I-Idir))/np.max(np.abs(Idir))
        print(f"{name}: rel err={err:.2e}   I[0:3]={I[:3]}  Idir[0:3]={Idir[:3]}")

if __name__=="__main__":
    main()