import numpy as np, time, sys, json
from global_bvp_g import solve_global
from dissipation_g import fields_T, gamma_D_ratio
from dissipation_psi import gamma_D_psi

def Kfit(X,T,lo,hi,order=4):
    m=(X>lo)&(X<hi); x=X[m]; t=T[m]
    A=np.column_stack([x**(-0.5),np.ones_like(x)]+[x**(k+0.5) for k in range(order)])
    c,*_=np.linalg.lstsq(A,t,rcond=None); return c

def diag(V,a,eps,Ny=320,N2=900,N3=2400,N1=200,L=16):
    r=solve_global(a,eps,V,L=L,N1=N1,N2=N2,N3=N3,Ny=Ny)
    X=r['X'];U0=r['U0'];T2=r['T2'];Y=r['Y'];U=r['U']
    it=np.argmin(np.abs(X+eps)); ue=U0[it+1]
    i0=np.argmin(np.abs(X))
    print('V',V,'a',a,'eps',eps,'Ny',Ny,'N2',N2,'N1',N1,'L',L)
    print('ue',ue,'2ue',2*ue,'Fcrit',2*ue-0.2,'T2 at i0',T2[i0],'X[i0]',X[i0])
    for hi in (0.03,0.05,0.1,0.2,0.4):
        c=Kfit(X,T2,0.003,hi,4)
        print(' Kfit hi %.3f C %+.3e D %+.5f E %+.3e'%(hi,c[0],c[1],c[2]))
    # T2 as function x
    for xx in (0,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,0.1,0.3):
        ii=np.argmin(np.abs(X-xx)); print('  T x %.4g X %.4g T %.6f'%(xx,X[ii],T2[ii]))
    gpsi=gamma_D_psi(r,V); gold=gamma_D_ratio(r,V);R=a*a/0.2
    print('gpsi %.8f gold %.8f R %.8f R-1-gpsi %.6f R-1-gold %.6f relpsi %.5f relold %.5f'%(gpsi,gold,R,R-1-gpsi,R-1-gold,(R-1-gpsi)/max(1,R),(R-1-gold)/max(1,R)))
    # residual at X=0 interface T2 jump
    print('T2 just left/right of 0')
    for ii in range(max(1,i0-3),min(len(X),i0+4)):
        print('   i %d X %+.6g T2 %.8f'%(ii,X[ii],T2[ii]))
    return r,ue

if __name__=='__main__':
    V=float(sys.argv[1]); a=float(sys.argv[2]); eps=float(sys.argv[3])
    kw={}
    for arg in sys.argv[4:]:
        k,v=arg.split('='); kw[k]=int(v) if k in ('Ny','N2','N3','N1') else float(v)
    diag(V,a,eps,**kw)
