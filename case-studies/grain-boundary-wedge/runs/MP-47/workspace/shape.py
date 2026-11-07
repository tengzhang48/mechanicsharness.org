import numpy as np
from scipy.linalg import null_space, eig
from explore_eig import clencurt, chebD, build

N=200
w0,xc=clencurt(N); x=(xc+1)/2.0; w=w0/2.0
D0,_=chebD(N); Dz=2*D0; D1=Dz; D2=Dz@Dz; D3=Dz@Dz@Dz
Aint=build(x,D1,D3,w,N)
C=np.vstack([D1[0,:],D2[0,:]]); B=null_space(C)
lam,vec=eig(Aint@B,B[1:N,:])
# find the converged mode near -106
j=np.argmin(np.abs(lam.real+106.0))
print("lam=",lam[j])
phi=np.real(B@vec[:,j])
print("phi(0)=%.3e phi(node1)=%.3e phi(1)=%.3e"%(phi[N],phi[N-1],phi[0]))
print("x (from 1 to 0) vs phi:")
for k in list(range(0,12))+list(range(N-8,N+1)):
    print("  x=%.6f  phi=%.6f"%(x[k],phi[k]/phi[0]))
# also normalized shape at x=0.5
idx=np.argmin(np.abs(x-0.5))
print("x[%d]=%.6f phi_rel=%.6f"%(idx,x[idx],phi[idx]/phi[0]))