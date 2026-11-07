import numpy as np, time
from global_bvp import solve_global
from dissipation import gamma_D_ratio

def state(a,eps,V,grid):
    r=solve_global(a,eps,V,**grid)
    X=r['X']; U0=r['U0']; it=np.argmin(np.abs(X+eps))
    ue=U0[it+1]
    g=gamma_D_ratio(r,V) if V!=0 else 0.0
    return ue,g

def extrap(a,eps,V,base,Nys=(200,400)):
    vals=[]
    for Ny in Nys:
        kw=dict(base); kw['Ny']=Ny
        ue,g=state(a,eps,V,kw)
        vals.append((1.0/Ny,ue,g))
    vals=np.array(vals)
    A=np.column_stack([np.ones(len(vals)),vals[:,0]])
    cu,*_=np.linalg.lstsq(A,vals[:,1],rcond=None)
    cg,*_=np.linalg.lstsq(A,vals[:,2],rcond=None)
    return cu[0],cg[0],vals

if __name__=='__main__':
    base=dict(N1=120,N2=360,N3=300)
    for V,a,eps in [(0,0.4472136,0.5),(0.1,0.5,0.5),(0.1,0.55,0.55)]:
        t0=time.time(); ue,g,vals=extrap(a,eps,V,base,Nys=(150,300,600))
        print('V %.3g a %.4f eps %.3f ue_inf %.6f g_inf %.6f  raw'%(V,a,eps,ue,g),vals[:,1:].T,'t %.1f'%(time.time()-t0))
