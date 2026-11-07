import numpy as np, time, sys
from global_bvp_g import solve_global
from scipy.optimize import least_squares

def Kfit(X,T,lo=0.01,hi=0.4,order=5):
    m=(X>lo)&(X<hi); x=X[m]; t=T[m]
    A=np.column_stack([x**(-0.5),np.ones_like(x)]+[x**(k+0.5) for k in range(order)])
    c,*_=np.linalg.lstsq(A,t,rcond=None); return c

def evalx(x,V,grid):
    a,eps=x
    r=solve_global(a,eps,V,**grid)
    X=r['X'];U0=r['U0'];it=np.argmin(np.abs(X+eps));ue=U0[it+1]
    c=Kfit(X,r['T2'])
    return ue,c[0],r

def res(x,V,grid):
    ue,C,r=evalx(x,V,grid)
    # scale K residual relative to needed state; C target 0. Use C*10 maybe
    return np.array([2*ue-0.2, C*10.0])

if __name__=='__main__':
    V=float(sys.argv[1]); x0=[float(sys.argv[2]),float(sys.argv[3])]
    grid=dict(L=16,N1=200,N2=900,N3=2400,Ny=320)
    t0=time.time()
    sol=least_squares(res,x0,args=(V,grid),bounds=([0.3,0.1],[1.5,2.5]),xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=100,x_scale=[0.03,0.05],verbose=1)
    print('SOL',V,sol.x,sol.fun,sol.nfev,'t',time.time()-t0)
    ue,C,r=evalx(sol.x,V,grid)
    print('DETAIL ue %.9f C %+.6e'%(ue,C))
