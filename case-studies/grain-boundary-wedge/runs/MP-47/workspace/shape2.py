import numpy as np
from scipy.linalg import null_space, eig
from explore_eig import clencurt, chebD
from op2 import build_Op2

N=400; ngl=700; drop='start'
w0,xc=clencurt(N); x=(xc+1)/2.0
D0,_=chebD(N); Dz=2*D0; D1=Dz; D2=Dz@Dz; D3=Dz@Dz@Dz
Op=build_Op2(x,N,ngl)
Aint=-Op@D3
eN=np.zeros(N+1); eN[N]=1.0
C=np.vstack([eN, D1[0,:], D2[0,:]])
B=null_space(C)
A=Aint@B; Mall=B[1:N,:]
keep=np.arange(1,N-1)
lam,vec=eig(A[keep,:], Mall[keep,:])
idx=np.argsort(lam.real)
for j in range(3):
    k=idx[j]
    print("lam=%.6f (imag %.2e)"%(lam[k].real, lam[k].imag))
c=vec[:,idx[0]].real
v=B@c
v=v/np.abs(v).max()
phi=v/v[0]                      # f(1)=1
node_half=N//2                  # x=0.5 for even N
print("x[half]=",x[node_half])
print("f(0.5)=",phi[node_half])
print("phi(1)=",phi[0]," phi(0)=",phi[N])
print("phi''(1)=%.3e phi'(1)=%.3e"%( (D2@v)[0]/np.abs(v).max(), (D1@v)[0]/np.abs(v).max() ))
print("phi'''(1)=%.3e"%( (D3@v)[0]/np.abs(v).max() ))
print("phi'''(0)=%.3e"%( (D3@v)[N]/np.abs(v).max() ))
# sample profile
for xv in [0.0,0.1,0.25,0.5,0.75,0.9,0.99,1.0]:
    i=int(np.argmin(np.abs(x-xv)))
    print("  x=%.2f  f=%.6f"%(x[i],phi[i]))