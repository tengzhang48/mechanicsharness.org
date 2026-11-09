"""Compare angular-form I_n vs Cartesian contour I_n term by term."""
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
    sepp_noD3=(se2pp-2.0*sep*sep)/(2.0*se)
    sepp_d3=se2pp_d3/(2.0*se)
    qpp_noD3=(n-1.0)*((n-2.0)*se**(n-3.0)*sep*sep+se**(n-2.0)*sepp_noD3)
    qpp_d3=(n-1.0)*se**(n-2.0)*sepp_d3
    RHS=(p+1.0)*((3.0/4.0)*q*(B-A)-Ur)+3.0*(qp*C+q*Cp)
    lhs_noD3=qpp_noD3*UmB+2.0*qp*UmB_p+q*(k*d2-Bpp)
    lhs_d3=qpp_d3*UmB+q
    d3p=((4.0/3.0)*p*RHS-lhs_noD3)/lhs_d3
    return [d1,d2,d3,d3p]

Y0=[1.0,0.0,-0.7050162232,0.0932567960]
s=solve_ivp(rhs,[0.0,np.pi],Y0,rtol=1e-13,atol=1e-15,method='Radau',dense_output=True)
sol=s.sol
tt=np.linspace(0,np.pi,40001)
Yv=sol(tt); ph,d1,d2,d3=Yv
sev=np.sqrt(0.75*(k*ph+d2-(k*(k-1.0)*ph))**2+3.0*(-(k-1.0)*d1)**2)
sc=1.0/sev.max()
def fields(thd):
    Y=sol(thd); ph,d1,d2,d3=Y
    A=(k*ph+d2)*sc; B=(k*(k-1.0)*ph)*sc; C=(-(k-1.0)*d1)*sc
    se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
    q=se**(n-1.0)
    return A,B,C,se,q
A,B,C,se,q=fields(tt)
err=0.75*q*(A-B); ett=0.75*q*(B-A); ert=1.5*q*C
Ur=0.75*q*(A-B)/p
Ut=(np.gradient(Ur,tt)-3.0*q*C)/(p+1.0)
def mir(x): return np.concatenate([x[::-1][:-1],x])
def miro(x): return np.concatenate([-x[::-1][:-1],x])  # odd extension
thf=np.concatenate([-tt[::-1][:-1],tt])
Af,Bf,Cf,errf,ettf,ertf,Urf=[mir(x) for x in (A,B,C,err,ett,ert,Ur)]
Utf=miro(Ut)   # u_theta is ODD in theta (Mode I symmetry)
Wt=(n/(n+1.0))*(Af*errf+Bf*ettf+2.0*Cf*ertf)
# angular J integrand
dUrf=np.gradient(Urf,thf); dUtf=np.gradient(Utf,thf)
# du/dx in polar->cart. du_i/dx = cos*u_{i,r} - sin/r * u_{i,th}
# u_r = C0 r^p Ur ; u_{r,r}=C0 p r^{p-1}Ur ; u_{r,th}=C0 r^p Ur'
# => u_{r,x}=C0 r^{p-1}(cos*p*Ur - sin*Ur')
# u_th component similarly
# radial comp of du/dx: cos*du_r/dr - sin/r*(du_r/dth - u_th)
#   u_r=C0 r^p Ur ; du_r/dr=C0 p r^{p-1} Ur ; du_r/dth=C0 r^p Ur'
#   => u_{r,x}=C0 r^{p-1}[ cos*p*Ur - sin*(Ur' - Ut) ]
# theta comp: cos*du_th/dr - sin/r*(du_th/dth + u_r)
#   => u_{th,x}=C0 r^{p-1}[ cos*p*Ut - sin*(Ut' + Ur) ]
ur_x=np.cos(thf)*p*Urf-np.sin(thf)*(dUrf-Utf)
uth_x=np.cos(thf)*p*Utf-np.sin(thf)*(dUtf+Urf)
# t_r=A, t_th=C ; t_i u_{i,x} = A*ur_x + C*uth_x
tt_int=Af*ur_x+Cf*uth_x
g = Wt*np.cos(thf)-tt_int
I_n=np.trapezoid(g,thf)
print("angular I_n =", I_n)
# Now integrand breakdown
print("int W cos =", np.trapezoid(Wt*np.cos(thf),thf))
print("int A*ur_x =", np.trapezoid(Af*ur_x,thf))
print("int C*uth_x =", np.trapezoid(Cf*uth_x,thf))