"""
CLEAN final computation with CORRECT consistent scaling.
Normalization: max over theta of se_tilde = 1.
Fields: sigma = sigma0 K r^{-p} sigtilde ; eps = alpha eps0 K^n r^{p-1} epstilde ;
        u = alpha eps0 K^n r^{p} utilde ; W = sigma0 alpha eps0 K^{n+1} r^{-1} Wtilde.
Mathematics fixed at: se_tilde^max = 1.
"""
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
tt=np.linspace(0,np.pi,200001)
Yv=sol(tt); ph,d1,d2,d3=Yv
Av=k*ph+d2; Bv=k*(k-1.0)*ph; Cv=-(k-1.0)*d1
sev=np.sqrt(0.75*(Av-Bv)**2+3.0*Cv**2)
sc=1.0/sev.max()   # scale factor making max se_tilde = 1
print("scale factor sc=%.10f (raw se_max=%.10f)" % (sc, sev.max()))

# scaled fields on [0,pi]
A=(k*ph+d2)*sc; B=(k*(k-1.0)*ph)*sc; C=(-(k-1.0)*d1)*sc
se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
q=se**(n-1.0)
err=0.75*q*(A-B); ett=0.75*q*(B-A); ert=1.5*q*C
# displacement
Ur=0.75*q*(A-B)/p
Ut=(np.gradient(Ur,tt)-3.0*q*C)/(p+1.0)
# mirror to [-pi,pi]
def mir(x): return np.concatenate([x[::-1][:-1],x])
thf=np.concatenate([-tt[::-1][:-1],tt])
Af,Bf,Cf,errf,ettf,ertf,Urf,Utf=[mir(x) for x in (A,B,C,err,ett,ert,Ur,Ut)]
dUrf=np.gradient(Urf,thf); dUtf=np.gradient(Utf,thf)
Wt=(n/(n+1.0))*(Af*errf+Bf*ettf+2.0*Cf*ertf)
dur_dx=np.cos(thf)*p*Urf-np.sin(thf)*dUrf+np.sin(thf)*Utf
dut_dx=np.cos(thf)*p*Utf-np.sin(thf)*dUtf-np.sin(thf)*Urf
tt_int=Af*dur_dx+Cf*dut_dx
I_n=np.trapezoid(Wt*np.cos(thf)-tt_int,thf)
Wbar=np.trapezoid(Wt,thf)
B0=B[0]; A0=A[0]
print("I_n=%.8f  Wbar=%.8f  B(0)=%.8f  A(0)=%.8f" % (I_n,Wbar,B0,A0))
print("H=%.8f  C=%.8f  p=%.8f" % (B0/I_n**p, Wbar/I_n, p))
# sanity: W tilde = se^{n+1}
print("max|Wtilde - se^{n+1}| =", np.max(np.abs(Wt-se**(n+1))))