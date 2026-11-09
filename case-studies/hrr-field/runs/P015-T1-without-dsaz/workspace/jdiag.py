"""Diagnose the J-integrand terms for the corrected solution."""
import numpy as np
from scipy.integrate import solve_ivp

n = 6.5
p = 1.0/(n+1.0)
k = 2.0-p


def rhs(theta, Y):
    phi, d1, d2, d3 = Y
    A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
    se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
    Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
    sep=(0.75*2.0*(A-B)*(Ap-Bp)+3.0*2.0*C*Cp)/(2.0*se)
    q=se**(n-1.0); qp=(n-1.0)*se**(n-2.0)*sep
    Ur=q*(A-B)/(2.0*p)
    Urp=(1.0/(2.0*p))*(qp*(A-B)+q*(Ap-Bp))
    Urpp=(p+1.0)*(q*(B-A)/2.0-Ur)+2.0*(qp*C+q*Cp)
    Bpp=k*(k-1.0)*d2; Cpp=-(k-1.0)*d3
    UmB=A-B; UmB_p=Ap-Bp; UmB_pp_noD3=k*d2-Bpp
    se2pp=0.75*(2.0*UmB_p*UmB_p+2.0*UmB*UmB_pp_noD3)+3.0*(2.0*Cp*Cp+2.0*C*Cpp)
    se2pp_d3=1.5*UmB
    sepp_noD3=(se2pp-2.0*sep*sep)/(2.0*se)
    sepp_d3=se2pp_d3/(2.0*se)
    qpp_noD3=(n-1.0)*((n-2.0)*se**(n-3.0)*sep*sep+se**(n-2.0)*sepp_noD3)
    qpp_d3=(n-1.0)*se**(n-2.0)*sepp_d3
    lhs_noD3=qpp_noD3*UmB+2.0*qp*UmB_p+q*(k*d2-Bpp)
    lhs_d3=qpp_d3*UmB+q
    d3p=(2.0*p*Urpp-lhs_noD3)/lhs_d3
    return [d1,d2,d3,d3p]


Y0=[1.0,0.0,-0.7050162232,0.0932567960]
s=solve_ivp(rhs,[0.0,np.pi],Y0,rtol=1e-12,atol=1e-14,method='Radau',dense_output=True,max_step=0.005)
th=s.t; Y=s.y.T
phi,d1,d2,d3=Y.T
A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
sep=(0.75*2.0*(A-B)*(Ap-Bp)+3.0*2.0*C*Cp)/(2.0*se)
q=se**(n-1.0); qp=(n-1.0)*se**(n-2.0)*sep
Ur=q*(A-B)/(2.0*p)
Urp=(1.0/(2.0*p))*(qp*(A-B)+q*(Ap-Bp))
Ut=(Urp-2.0*q*C)/(p+1.0)
sc=1.0/se.max()
A,B,C,se,Ur,Urp,Ut=[x*sc for x in (A,B,C,se,Ur,Urp,Ut)]
def mir(x): return np.concatenate([x[::-1][:-1],x])
th_f=np.concatenate([-th[::-1][:-1],th])
Af,Bf,Cf=mir(A),mir(B),mir(C)
Urf,Utf,Urpf=mir(Ur),mir(Ut),mir(Urp)
dUtf=np.gradient(Utf,th_f)
err=p*Urf; ett=Urf+dUtf; ert=0.5*(Urpf-p*Utf-Utf)
Wt=(n/(n+1.0))*(Af*err+Bf*ett+2.0*Cf*ert)
dur_dx=np.cos(th_f)*p*Urf-np.sin(th_f)*Urpf+np.sin(th_f)*Utf
dut_dx=np.cos(th_f)*p*Utf-np.sin(th_f)*dUtf-np.sin(th_f)*Urf
tt=Af*dur_dx+Cf*dut_dx
print("integrand W cos th - t.u_x : sample (deg)")
for i in range(0, len(th_f), len(th_f)//12):
    print("  th=%6.1f  Wcos=%.5f  tt=%.5f  diff=%.5f" %
          (np.degrees(th_f[i]), Wt[i]*np.cos(th_f[i]), tt[i],
           Wt[i]*np.cos(th_f[i])-tt[i]))
print("I_n =", np.trapezoid(Wt*np.cos(th_f)-tt, th_f))
print("int Wt cos =", np.trapezoid(Wt*np.cos(th_f), th_f))
print("int tt =", np.trapezoid(tt, th_f))
print("Wbar =", np.trapezoid(Wt, th_f))