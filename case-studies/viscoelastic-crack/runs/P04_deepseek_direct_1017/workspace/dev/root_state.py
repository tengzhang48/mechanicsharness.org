import numpy as np, time, json, sys
from global_bvp import solve_global
from dissipation import gamma_D_ratio

DELTA=0.2

def eval_state(a,eps,V,grid):
    r=solve_global(a,eps,V,**grid)
    X=r['X']; U0=r['U0']
    it=np.argmin(np.abs(X+eps))
    ue=U0[it+1]
    g=gamma_D_ratio(r,V) if V!=0 else 0.0
    return r,ue,g

def residual(x,V,grid,verbose=False):
    a,eps=x
    r,ue,g=eval_state(a,eps,V,grid)
    F1=2*ue-DELTA
    F2=a*a/DELTA - 2*ue/DELTA - g
    if verbose: print('  a %.6f eps %.6f ue %.6f g %.5f F1 %.5f F2 %.5f'%(a,eps,ue,g,F1,F2))
    return np.array([F1,F2])

if __name__=='__main__':
    from scipy.optimize import root
    grids=[dict(N1=100,N2=300,N3=240,Ny=200)]
    for V in [float(v) for v in (sys.argv[1:] or ['0'])]:
        for grid in grids:
            t0=time.time()
            sol=root(residual,[0.45,0.5],args=(V,grid),method='hybr',options=dict(xtol=1e-8))
            print('V',V,'grid',grid,'success',sol.success,'x',sol.x,'fun',sol.fun,'time %.1f'%(time.time()-t0))
            residual(sol.x,V,grid,verbose=True)
