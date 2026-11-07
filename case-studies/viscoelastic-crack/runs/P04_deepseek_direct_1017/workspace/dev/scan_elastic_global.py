import numpy as np, time
from global_bvp import solve_global

def Kfit(X,T2):
    # fit T2 = C/sqrt(X)+D+E sqrt(X) on near-tip points
    m=(X>0)&(X<0.15)
    x=X[m]; t=T2[m]
    if len(x)<8: return np.nan
    A=np.column_stack([1/np.sqrt(x),np.ones_like(x),np.sqrt(x)])
    coef,*_=np.linalg.lstsq(A,t,rcond=None)
    return coef[0]

if __name__=='__main__':
    a=np.sqrt(0.2)
    for eps in [0.45,0.48,0.50,0.515,0.53,0.56,0.60]:
        r=solve_global(a,eps,0.0)
        X=r['X']; U0=r['U0']; T2=r['T2']
        ie=np.argmin(np.abs(X+eps)); i0=np.argmin(np.abs(X))
        print('eps %.3f U0(-eps) %.5f T2(0+) %.4f K %.4f'%(eps,U0[ie],T2[i0],Kfit(X,T2)))
