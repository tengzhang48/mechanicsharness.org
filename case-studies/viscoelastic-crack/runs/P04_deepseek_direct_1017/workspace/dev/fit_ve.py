import numpy as np, time
from global_bvp_g import solve_global
from dissipation_g import gamma_D_ratio
D=0.2

def res(x,V,kw):
    a,eps=x
    if a<=0 or eps<=0: return np.array([1e3,1e3])
    r=solve_global(a,eps,V,**kw)
    X=r['X']; U0=r['U0']; it=np.argmin(np.abs(X+eps)); ue=U0[it+1]
    g=gamma_D_ratio(r,V) if V else 0.0
    F1=2*ue-D
    F2=a*a/D-2*ue/D-g
    return np.array([F1,F2])

if __name__=='__main__':
    from scipy.optimize import least_squares
    kw=dict(N1=100,N2=300,N3=240,Ny=120)
    V=0.1
    for x0 in [[0.45,0.5],[0.4,0.3],[0.5,0.7],[0.42,0.2],[0.6,1.0],[0.38,0.15],[0.55,1.5]]:
        t0=time.time()
        sol=least_squares(res,x0,args=(V,kw),bounds=([0.01,0.01],[3,5]),xtol=1e-10,ftol=1e-10,gtol=1e-10,max_nfev=200)
        print('x0',x0,'->',sol.x,'res',sol.fun,'cost',sol.cost,'nfev',sol.nfev,'t %.1f'%(time.time()-t0))
