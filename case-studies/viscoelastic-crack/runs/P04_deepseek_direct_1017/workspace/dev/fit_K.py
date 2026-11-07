import numpy as np, time
from global_bvp_g import solve_global

def Kcoef(X,T,lo=0.05,hi=0.4):
    m=(X>lo)&(X<hi); x=X[m]; t=T[m]
    A=np.column_stack([1/np.sqrt(x),np.ones_like(x),np.sqrt(x),x])
    c,*_=np.linalg.lstsq(A,t,rcond=None); return c[0]

def resK(x,V,kw):
    a,eps=x
    r=solve_global(a,eps,V,**kw)
    X=r['X'];U0=r['U0'];it=np.argmin(np.abs(X+eps));ue=U0[it+1]
    C=Kcoef(X,r['T2'])
    return np.array([2*ue-0.2, C/0.01])

if __name__=='__main__':
    import sys
    from scipy.optimize import least_squares
    kw=dict(L=16,N1=120,N2=360,N3=600,Ny=160)
    V=float(sys.argv[1]); x0=[float(sys.argv[2]),float(sys.argv[3])]
    t0=time.time()
    sol=least_squares(resK,x0,args=(V,kw),bounds=([0.3,0.1],[1.5,2.5]),xtol=1e-10,ftol=1e-10,gtol=1e-10,max_nfev=80,x_scale=[0.05,0.1])
    print('V',V,'x0',x0,'->',sol.x,'fun',sol.fun,'nfev',sol.nfev,'t %.1f'%(time.time()-t0))
    a,eps=sol.x
    r=solve_global(a,eps,V,**kw)
    X=r['X'];U0=r['U0'];it=np.argmin(np.abs(X+eps))
    print('ue',U0[it+1],'C',Kcoef(X,r['T2']))
