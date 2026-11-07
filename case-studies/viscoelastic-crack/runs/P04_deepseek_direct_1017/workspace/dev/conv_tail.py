import numpy as np, time
from global_bvp import solve_global
from tipfit2 import tipfit

def meas(eps,a,V=0.0,**kw):
    r=solve_global(a,eps,V,**kw)
    X=r['X']; T=r['T2']; U0=r['U0']
    it=np.argmin(np.abs(X+eps))
    ue=U0[it+1]
    D=tipfit(X,T,0.03)[1]
    C=tipfit(X,T,0.03)[0]
    return ue,D,C

if __name__=='__main__':
    eps=0.523; a=np.sqrt(0.2)
    for kw in [dict(N1=120,N2=360,N3=300,Ny=200),
               dict(N1=160,N2=480,N3=400,Ny=300),
               dict(N1=200,N2=600,N3=500,Ny=400),
               dict(N1=240,N2=720,N3=600,Ny=500)]:
        t0=time.time(); ue,D,C=meas(eps,a,**kw)
        print(kw,'U0 %.6f D %.6f C %.6f  t %.1f'%(ue,D,C,time.time()-t0))
