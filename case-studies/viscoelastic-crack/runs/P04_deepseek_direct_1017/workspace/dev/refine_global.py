import numpy as np
from global_bvp import solve_global
from scan_elastic_global import Kfit

def meas(eps,V=0.0,a=np.sqrt(0.2),**kw):
    r=solve_global(a,eps,V,**kw)
    X=r['X']; U0=r['U0']; T2=r['T2']
    ie=np.argmin(np.abs(X+eps)); i0=np.argmin(np.abs(X))
    return U0[ie],T2[i0],Kfit(X,T2),X,T2,U0

if __name__=='__main__':
    for kw in [dict(N1=80,N2=240,N3=120,Ny=150),
               dict(N1=120,N2=360,N3=180,Ny=200),
               dict(N1=160,N2=480,N3=240,Ny=250),
               dict(N1=200,N2=600,N3=300,Ny=300)]:
        ue,t0,k,_,_,_=meas(0.5,**kw)
        print(kw,'U0 %.6f T2(0+) %.6f K %.5f'%(ue,t0,k))
    print('--- eps scan fine ---')
    for eps in np.arange(0.49,0.551,0.01):
        ue,t0,k,_,_,_=meas(eps,N2=360,N3=180,Ny=200,N1=120)
        print('eps %.3f U0 %.6f T2(0+) %.6f K %.5f'%(eps,ue,t0,k))
