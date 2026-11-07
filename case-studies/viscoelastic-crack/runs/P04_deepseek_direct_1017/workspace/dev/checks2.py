import numpy as np
from global_bvp_g import solve_global
from dissipation_g import fields_T
from dissipation_psi import gamma_D_psi
from scipy.interpolate import interp1d

def Kfit(X,T,lo=.01,hi=.4,order=5):
    m=(X>lo)&(X<hi);x=X[m];t=T[m]
    A=np.column_stack([x**(-.5),np.ones_like(x)]+[x**(k+.5) for k in range(order)])
    c,*_=np.linalg.lstsq(A,t,rcond=None); return c

def run(V,a,eps):
    grid=dict(L=16,N1=200,N2=900,N3=2400,Ny=320)
    r=solve_global(a,eps,V,**grid); X=r['X'];Y=r['Y'];U=r['U'];Nx=len(X)
    E1,E2,T1,T2=fields_T(r,V)
    def dX(F):
        out=np.zeros_like(F)
        for k in range(F.shape[1]):
            i=k+1
            if i-1>=1 and i+1<=Nx-2:
                hm=X[i]-X[i-1];hp=X[i+1]-X[i]
                out[:,k]=-hp/(hm*(hm+hp))*F[:,k-1]+(hp-hm)/(hm*hp)*F[:,k]+hm/(hp*(hm+hp))*F[:,k+1]
            elif i-1<1: out[:,k]=(F[:,k+1]-F[:,k])/(X[i+1]-X[i])
            else: out[:,k]=(F[:,k]-F[:,k-1])/(X[i]-X[i-1])
        return out
    dT1=dX(T1)
    dT2y=np.zeros_like(T2)
    for j in range(T2.shape[0]):
        jj=j+1
        if jj==1: dT2y[j,:]=(T2[j+1,:]-T2[j,:])/(Y[2]-Y[1])
        elif jj==len(Y)-2: dT2y[j,:]=(T2[j,:]-T2[j-1,:])/(Y[-2]-Y[-3])
        else:
            hm=Y[jj]-Y[jj-1];hp=Y[jj+1]-Y[jj]
            dT2y[j,:]=-hp/(hm*(hm+hp))*T2[j-1,:]+(hp-hm)/(hm*hp)*T2[j,:]+hm/(hp*(hm+hp))*T2[j+1,:]
    res=dT1+dT2y; denom=max(1,np.sqrt(np.mean(dT1**2+dT2y**2)))
    Yi=Y[1:-1]
    for lo in [0.0,0.01,0.05,0.3,0.6]:
        m=Yi>=lo
        val=np.sqrt(np.mean(res[m,:]**2))/denom
        print('V',V,'eq Y>',lo,'rel',val)
    c=Kfit(X,r['T2']); print('V',V,'Kfit C %.4e D %.6f'%(c[0],c[1]))
    i0=np.argmin(abs(X)); it=np.argmin(abs(X+eps))
    print('V',V,'ue',r['U0'][it+1],'crit',2*r['U0'][it+1]-.2,'T2(0+)',r['T2'][i0])
    print('V',V,'energy',a*a/.2-1-gamma_D_psi(r,V),'R',a*a/.2,'g',gamma_D_psi(r,V))
for V,a,eps in [(0.001,.451,.5489),(0.1,.5175,.755)]: run(V,a,eps)
