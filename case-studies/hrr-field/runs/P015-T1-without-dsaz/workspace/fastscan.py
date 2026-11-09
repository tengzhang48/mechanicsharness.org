"""Fast scan of the (phi''(0),phi'''(0)) plane for zeros of (A(pi),C(pi))."""
import numpy as np
from scipy.integrate import solve_ivp

n=6.5; p=1.0/(n+1.0); k=2.0-p
def rhs(theta,Y):
    phi,d1,d2,d3=Y
    A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
    se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
    Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
    sep=(0.75*2.0*(A-B)*(Ap-Bp)+3.0*2.0*C*Cp)/(2.0*se)
    q=se**(n-1.0); qp=(n-1.0)*se**(n-2.0)*sep
    Ur=(3.0/4.0)*q*(A-B)/p
    Bpp=k*(k-1.0)*d2; Cpp=-(k-1.0)*d3
    UmB=A-B; UmB_p=Ap-Bp; UmB_pp_noD3=k*d2-Bpp
    se2pp=0.75*(2.0*UmB_p*UmB_p+2.0*UmB*UmB_pp_noD3)+3.0*(2.0*Cp*Cp+2.0*C*Cpp)
    se2pp_d3=1.5*UmB
    sepp_noD3=(se2pp-2.0*sep*sep)/(2.0*se); sepp_d3=se2pp_d3/(2.0*se)
    qpp_noD3=(n-1.0)*((n-2.0)*se**(n-3.0)*sep*sep+se**(n-2.0)*sepp_noD3)
    qpp_d3=(n-1.0)*se**(n-2.0)*sepp_d3
    RHS=(p+1.0)*((3.0/4.0)*q*(B-A)-Ur)+3.0*(qp*C+q*Cp)
    lhs_noD3=qpp_noD3*UmB+2.0*qp*UmB_p+q*(k*d2-Bpp); lhs_d3=qpp_d3*UmB+q
    d3p=((4.0/3.0)*p*RHS-lhs_noD3)/lhs_d3
    return [d1,d2,d3,d3p]

def resid(u):
    s=solve_ivp(rhs,[0.0,np.pi],[1.0,0.0,u[0],u[1]],rtol=1e-7,atol=1e-9,
                method='RK45',dense_output=False)
    if not s.success: return np.array([np.nan,np.nan])
    phi,d1,d2,d3=s.y[:,-1]
    return np.array([k*phi+d2, -(k-1.0)*d1])

a=np.linspace(-1.6,0.5,80); b=np.linspace(-0.5,0.6,80)
F=np.zeros((len(a),len(b),2))
for i,ai in enumerate(a):
    for j,bj in enumerate(b):
        F[i,j]=resid([ai,bj])
nA=F[:,:,0]; nC=F[:,:,1]
# find sign-change cells in both components simultaneously
cands=[]
for i in range(len(a)-1):
    for j in range(len(b)-1):
        def signchange(M):
            v=M[i:i+2,j:j+2].ravel()
            return np.any(v>0) and np.any(v<0)
        if signchange(nA) and signchange(nC):
            cands.append((0.5*(a[i]+a[i+1]),0.5*(b[j]+b[j+1])))
print("candidate cells:",len(cands))
for c in cands[:40]:
    r=resid(c)
    print("  guess (%.3f,%.3f) -> resid (%.4f,%.4f)"%(c[0],c[1],r[0],r[1]))
np.savez('scangrid.npz',a=a,b=b,nA=nA,nC=nC)