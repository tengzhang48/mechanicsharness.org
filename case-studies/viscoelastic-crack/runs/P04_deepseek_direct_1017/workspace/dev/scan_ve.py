import numpy as np, time
from global_bvp import solve_global
from dissipation import gamma_D_ratio
D=0.2

def F(a,eps,V,kw):
    r=solve_global(a,eps,V,**kw)
    X=r['X']; U0=r['U0']; it=np.argmin(np.abs(X+eps))
    ue=U0[it+1]
    g=gamma_D_ratio(r,V) if V!=0 else 0.0
    F1=2*ue-D
    F2=a*a/D-2*ue/D-g
    return F1,F2,ue,g

if __name__=='__main__':
    kw=dict(N1=100,N2=300,N3=240,Ny=200)
    for V in [0.1]:
        print('V',V)
        for a in [0.45,0.5,0.55,0.6,0.65,0.7,0.8]:
            row=[]
            for eps in [0.4,0.5,0.6,0.7,0.8,1.0]:
                F1,F2,ue,g=F(a,eps,V,kw)
                row.append('%.3f/%.2f'%(ue,F2))
            print('a %.2f'%a, ' '.join(row))
