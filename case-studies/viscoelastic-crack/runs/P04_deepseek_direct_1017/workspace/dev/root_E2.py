import numpy as np, time, sys
from global_bvp_g import solve_global
from dissipation_g import gamma_D_ratio
from dissipation_psi import gamma_D_psi
from scipy.optimize import least_squares
D=.2

def evalx(x,V,grid):
    a,eps=x
    r=solve_global(a,eps,V,**grid)
    X=r['X'];U0=r['U0'];it=np.argmin(np.abs(X+eps));ue=U0[it+1]
    g=gamma_D_ratio(r,V) if V else 0.0
    gp=gamma_D_psi(r,V) if V else 0.0
    return ue,g,gp,r

def res(x,V,grid):
    a,eps=x
    ue,g,gp,r=evalx(x,V,grid)
    R=a*a/D
    return np.array([2*ue-D, (R-1-g)*5.0])

if __name__=='__main__':
    V=float(sys.argv[1]); x0=[float(sys.argv[2]),float(sys.argv[3])]
    grid=dict(L=16,N1=200,N2=900,N3=2400,Ny=320)
    t0=time.time()
    sol=least_squares(res,x0,args=(V,grid),bounds=([0.3,0.1],[1.5,2.5]),xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=100,x_scale=[0.03,0.05],verbose=1)
    print('SOL',V,sol.x,sol.fun,sol.nfev,'t',time.time()-t0)
    ue,g,gp,r=evalx(sol.x,V,grid);a,eps=sol.x;R=a*a/D
    print('DETAIL ue %.9f gold %.8f gpsi %.8f R %.8f R-1-gold %.5f rel %.5f R-1-gpsi %.5f rel %.5f'%(ue,g,gp,R,R-1-g,(R-1-g)/max(1,R),R-1-gp,(R-1-gp)/max(1,R)))
