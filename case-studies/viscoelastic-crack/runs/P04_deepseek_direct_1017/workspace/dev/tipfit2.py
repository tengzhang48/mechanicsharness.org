import numpy as np
from global_bvp import solve_global

def tipfit(X,T2,xmax=0.03):
    m=(X>0)&(X<xmax)
    x=X[m]; t=T2[m]
    A=np.column_stack([1/np.sqrt(x),np.ones_like(x),np.sqrt(x)])
    c,*_=np.linalg.lstsq(A,t,rcond=None)
    return c

def measure(a,eps,V,**kw):
    r=solve_global(a,eps,V,**kw)
    X=r['X']; T=r['T2']; U0=r['U0']
    it=np.argmin(np.abs(X+eps))
    # opening just inside cohesive zone
    ue=U0[it+1]
    C,D,E=tipfit(X,T,0.03)
    i0=np.argmin(np.abs(X))
    return ue,D,C,T[i0],r

if __name__=='__main__':
    a=np.sqrt(0.2)
    for eps in np.arange(0.48,0.556,0.005):
        ue,D,C,t0,r=measure(a,eps,0.0,N1=120,N2=360,N3=300,Ny=250)
        print('eps %.3f U0in %.6f D %.5f C %.5f T2(0+) %.5f'%(eps,ue,D,C,t0))
