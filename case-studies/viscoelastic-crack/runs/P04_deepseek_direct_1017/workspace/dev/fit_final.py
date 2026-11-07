import numpy as np, time
from global_bvp_g import solve_global
from dissipation_g import gamma_D_ratio
D=0.2

def ofit_c0(X,U0,lo=0.15,hi=0.002,order=2):
    m=(X<0)&(np.abs(X)>hi)&(np.abs(X)<lo)
    x=np.abs(X[m]); u=U0[m]
    A=np.column_stack([x**(0.5+k) for k in range(order+1)])
    c,*_=np.linalg.lstsq(A,u,rcond=None)
    return c[0]

def evalx(x,V,grid,full=False):
    a,eps=x
    r=solve_global(a,eps,V,**grid)
    X=r['X'];U0=r['U0'];it=np.argmin(np.abs(X+eps))
    ue=U0[it+1]; c0=ofit_c0(X,U0)
    if full:
        g=gamma_D_ratio(r,V) if V else 0.0
        return ue,c0,g,r
    return ue,c0

def res(x,V,grid):
    ue,c0=evalx(x,V,grid)
    return np.array([2*ue-D, 10*c0])

if __name__=='__main__':
    import sys
    from scipy.optimize import least_squares
    grid=dict(L=16,N1=200,N2=600,N3=1200,Ny=240)
    V=float(sys.argv[1]) if len(sys.argv)>1 else 0.0
    x0=[0.4472,0.541] if V==0 else ([0.4475,0.55] if V<0.01 else [0.5,0.6])
    t0=time.time()
    ue,c0,g,r=evalx(x0,V,grid,full=True)
    print('start',x0,'ue %.6f c0 %+.6f g %.5f F1 %.5f'%(ue,c0,g,2*ue-D))
    sol=least_squares(res,x0,args=(V,grid),bounds=([0.3,0.1],[1.5,2.5]),xtol=1e-10,ftol=1e-10,gtol=1e-10,max_nfev=100,x_scale=[0.05,0.1])
    print('sol',sol.x,'fun',sol.fun,'cost',sol.cost,'nfev',sol.nfev,'time %.1f'%(time.time()-t0))
    ue,c0,g,r=evalx(sol.x,V,grid,full=True)
    print('final ue %.6f c0 %+.6f g %.5f F1 %.5f R %.5f R-1-g %.5f'%(ue,c0,g,2*ue-D,sol.x[0]**2/D,sol.x[0]**2/D-1-g))
