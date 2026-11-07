import numpy as np, sys, time, json
from global_bvp_g import solve_global
from dissipation_g import fields_T
V=float(sys.argv[1]); a=float(sys.argv[2]); eps=float(sys.argv[3])
kw=dict(L=16,N1=200,N2=900,N3=2400,Ny=320)
r=solve_global(a,eps,V,**kw)
E1,E2,T1,T2=fields_T(r,V); X=r['X'];Y=r['Y'];Xi=X[1:-1];Yi=Y[1:-1]
Q2=(T1-E1)**2+(T2-E2)**2
# integral y for each x
fy=np.trapezoid(Q2,Yi,axis=0)
# cumulative x
Fx=np.concatenate([[0],np.cumsum((fy[1:]+fy[:-1])/2*np.diff(Xi))])
I=np.trapezoid(fy,Xi)
delta=.2; lam=.1
gold=2*I/(np.pi*V*(1-lam)*delta)
print('gold',gold,'I',I,'expected I',gold*np.pi*V*(1-lam)*delta/2)
print('tail left last (x<-10)',np.trapezoid(fy[Xi<-10],Xi[Xi<-10]),'frac',np.trapezoid(fy[Xi<-10],Xi[Xi<-10])/I)
print('tail right x>10',np.trapezoid(fy[Xi>10],Xi[Xi>10]),'frac',np.trapezoid(fy[Xi>10],Xi[Xi>10])/I)
print('fraction x<-4',np.trapezoid(fy[Xi<-4],Xi[Xi<-4])/I,'x>4',np.trapezoid(fy[Xi>4],Xi[Xi>4])/I)
for xx in (-16,-12,-8,-6,-4,-3,-2,-1,-.75,0,1,2,4,8,12,16):
    i=np.argmin(abs(Xi-xx)); print('x %+.3f fy %.6e Q2 max %.3e'%(Xi[i],fy[i],Q2[:,i].max()))
# Y layers contribution
gy=np.array([np.trapezoid(Q2[j,:],Xi) for j in range(len(Yi))])
print('Y frac first 5',gy[:5]/gy.sum(),'last5',gy[-5:]/gy.sum())
# integral near interface first dy contribution
print('first Y row frac of full trapezoid', (gy[0]*(Yi[1]-Yi[0])/2)/I)
np.savez('/workspace/dev/sol_%g_%g_%g.npz'%(V,a,eps),X=X,Y=Y,U=r['U'],E2i=r['E2'],T2i=r['T2'],U0=r['U0'])
