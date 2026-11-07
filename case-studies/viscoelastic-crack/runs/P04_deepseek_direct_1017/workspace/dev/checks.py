import numpy as np
from global_bvp_g import solve_global
from dissipation_g import fields_T
from dissipation_psi import gamma_D_psi

def checks(V,a,eps,grid):
    r=solve_global(a,eps,V,**grid)
    X=r['X'];Y=r['Y'];U=r['U'];dY=r['dY']
    E1,E2,T1,T2=fields_T(r,V)
    Xi=X[1:-1];Yi=Y[1:-1]
    it=np.argmin(np.abs(X+eps)); i0=np.argmin(np.abs(X))
    ue=r['U0'][it+1]
    # interface conditions
    mcoh=(X>=-eps)&(X<0); mcrack=(X<-eps)
    res_coh=np.max(np.abs(r['T2'][mcoh]-1.0))
    res_crack=np.max(np.abs(r['T2'][mcrack]))
    # material law residual in bulk
    def dX(F):
        out=np.zeros_like(F)
        for k in range(F.shape[1]):
            i=k+1
            if i-1>=1 and i+1<=len(X)-2:
                hm=X[i]-X[i-1];hp=X[i+1]-X[i]
                out[:,k]=-hp/(hm*(hm+hp))*F[:,k-1]+(hp-hm)/(hm*hp)*F[:,k]+hm/(hp*(hm+hp))*F[:,k+1]
            elif i-1<1:
                out[:,k]=(F[:,k+1]-F[:,k])/(X[i+1]-X[i])
            else:
                out[:,k]=(F[:,k]-F[:,k-1])/(X[i]-X[i-1])
        return out
    r1=T1-0.1*V*dX(T1)-(E1-V*dX(E1))
    r2=T2-0.1*V*dX(T2)-(E2-V*dX(E2))
    mat=np.sqrt(np.mean(r1**2+r2**2))/max(1,np.sqrt(np.mean(T1**2+T2**2+E1**2+E2**2)))
    # equilibrium residual tau1,x + tau2,y
    dT1=dX(T1)
    dT2y=np.zeros_like(T2)
    for j in range(T2.shape[0]):
        jj=j+1
        if jj==1: dT2y[j,:]=(T2[j+1,:]-T2[j,:])/(Y[2]-Y[1])
        elif jj==len(Y)-2: dT2y[j,:]=(T2[j,:]-T2[j-1,:])/(Y[-2]-Y[-3])
        else:
            hm=Y[jj]-Y[jj-1];hp=Y[jj+1]-Y[jj]
            dT2y[j,:]=-hp/(hm*(hm+hp))*T2[j-1,:]+(hp-hm)/(hm*hp)*T2[j,:]+hm/(hp*(hm+hp))*T2[j+1,:]
    eq=np.sqrt(np.mean((dT1+dT2y)**2))/max(1,np.sqrt(np.mean(dT1**2+dT2y**2)))
    g=gamma_D_psi(r,V) if V else 0.0
    R=a*a/0.2
    return dict(ue=ue,crit=ue-0.1,res_coh=res_coh,res_crack=res_crack,mat=mat,eq=eq,g=g,R=R,
                resid=R-1-g,rel=(R-1-g)/max(1,R),far_left=r['T2'][1],far_right=r['T2'][-2],
                Utop=np.max(np.abs(U[-1,:]-a)),Ulft=np.max(np.abs(U[:,0]-a)))

if __name__=='__main__':
    grid=dict(L=16,N1=200,N2=900,N3=2400,Ny=320)
    for V,a,eps in [(0.001,0.4510,0.5489),(0.1,0.5175,0.755)]:
        d=checks(V,a,eps,grid)
        print('V',V)
        for k,v in d.items(): print('  %-10s %+.3e'%(k,v))
