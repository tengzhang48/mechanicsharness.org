"""
Independent solve via solve_bvp with the corrected ODE and BCs:
  state [phi, phi', phi'', phi'''] on [0,pi]
  BCs: phi(0)=1, phi'(0)=0, A(pi)=0, C(pi)=0.
Then recompute I_n, Wbar, H, C.
"""
import numpy as np
from scipy.integrate import solve_bvp

n=6.5; p=1.0/(n+1.0); k=2.0-p

def ode(theta,Y):
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
    sepp_noD3=(se2pp-2.0*sep*sep)/(2.0*se)
    sepp_d3=se2pp_d3/(2.0*se)
    qpp_noD3=(n-1.0)*((n-2.0)*se**(n-3.0)*sep*sep+se**(n-2.0)*sepp_noD3)
    qpp_d3=(n-1.0)*se**(n-2.0)*sepp_d3
    RHS=(p+1.0)*((3.0/4.0)*q*(B-A)-Ur)+3.0*(qp*C+q*Cp)
    lhs_noD3=qpp_noD3*UmB+2.0*qp*UmB_p+q*(k*d2-Bpp)
    lhs_d3=qpp_d3*UmB+q
    d3p=((4.0/3.0)*p*RHS-lhs_noD3)/lhs_d3
    return np.vstack((d1,d2,d3,d3p))

def bc(Ya,Yb):
    # phi(0)=1, phi'(0)=0, A(pi)=0, C(pi)=0
    A_pi=k*Yb[0]+Yb[2]
    C_pi=-(k-1.0)*Yb[1]
    return np.array([Ya[0]-1.0, Ya[1], A_pi, C_pi])

x=np.linspace(0,np.pi,200)
y0=np.vstack((np.cos(x)*1.0, -np.sin(x), -np.cos(x), np.sin(x)))
res=solve_bvp(ode,bc,x,y0,tol=1e-11,max_nodes=20000)
print("solve_bvp success:", res.success, " nodes:", len(res.x))
Y=res.sol(res.x)
phi,d1,d2,d3=Y
A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
print("phi''(0)=%.7f phi'''(0)=%.7f" % (d2[0], d3[0]))
print("A(pi)=%.2e C(pi)=%.2e" % (A[-1], C[-1]))
se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
sc=1.0/se.max()
A=A*sc; B=B*sc; C=C*sc; se=se*sc
q=se**(n-1.0); Wt=(n/(n+1.0))*se**(n+1)
Ur=0.75*q*(A-B)/p; dUr=np.gradient(Ur,res.x); Ut=(dUr-3.0*q*C)/(p+1.0)
ur_x=np.cos(res.x)*p*Ur-np.sin(res.x)*(dUr-Ut)
g=Wt*np.cos(res.x)-A*ur_x
I_n=2.0*np.trapezoid(g,res.x)
Wbar=2.0*np.trapezoid(Wt,res.x)
print("I_n=%.7f Wbar=%.7f B(0)=%.7f" % (I_n,Wbar,B[0]))
print("H=%.7f C=%.7f" % (B[0]/I_n**p, Wbar/I_n))