import numpy as np
from global_bvp import solve_global

def tipfit(X,T2,xmax=0.03):
    m=(X>0)&(X<xmax)
    x=X[m]; t=T2[m]
    A=np.column_stack([1/np.sqrt(x),np.ones_like(x),np.sqrt(x)])
    c,*_=np.linalg.lstsq(A,t,rcond=None)
    return c

if __name__=='__main__':
    for eps in np.arange(0.47,0.561,0.01):
        r=solve_global(np.sqrt(0.2),eps,0.0,N1=120,N2=360,N3=300,Ny=200)
        X=r['X']; T=r['T2']; U0=r['U0']
        ie=np.argmin(np.abs(X+eps))
        C,D,E=tipfit(X,T)
        print('eps %.3f U0 %.5f  C %.5f D %.5f E %.5f'%(eps,U0[ie],C,D,E))
