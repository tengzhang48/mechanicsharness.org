"""Test teammate's (III) sign: U_r' - (p+1)U_t = 3qC, i.e. U_t=(U_r'-3qC)/(p+1).
With (I) pUr=(3/4)q(A-B), (II) Ur+Ut'=(3/4)q(B-A).
Derive d4. Check if an EVEN solution (phi'''(0)=0) then exists."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
n=6.5;p=1.0/(n+1);k=2-p

def rhs(theta,Y):
    phi,d1,d2,d3=Y
    A=k*phi+d2;B=k*(k-1)*phi;C=-(k-1)*d1
    se=np.sqrt(0.75*(A-B)**2+3*C**2);se=max(se,1e-30)
    Ap=k*d1+d3;Bp=k*(k-1)*d1;Cp=-(k-1)*d2
    UmB=A-B;UmB_p=Ap-Bp
    sep=(0.75*UmB*UmB_p+3*C*Cp)/se
    q=se**(n-1);qp=(n-1)*se**(n-2)*sep
    Ur=0.75*q*UmB/p
    # (III)-team: Ur'-(p+1)Ut=3qC => Ut=(Ur'-3qC)/(p+1)
    # (II): Ut'=(3/4)q(B-A)-Ur
    # diff(III): -(p+1)Ut'=3(qp C+q Cp)-Ur'' => Ut'=[Ur''-3(qpC+qCp)]/(p+1)
    # equate: [Ur''-3(qpC+qCp)]/(p+1) = (3/4)q(B-A)-Ur
    Urpp=3*(qp*C+q*Cp)+(p+1)*(0.75*q*(B-A)-Ur)
    Bpp=k*(k-1)*d2;Cpp=-(k-1)*d3
    S0=1.5*(UmB_p**2+UmB*(k*p*d2-Bpp))+6*(Cp**2+C*Cpp)
    sepp0=(S0-2*sep**2)/(2*se);seppd4=1.5*UmB/(2*se)
    qpp0=(n-1)*(n-2)*se**(n-3)*sep**2+(n-1)*se**(n-2)*sepp0
    qppd4=(n-1)*se**(n-2)*seppd4
    d4=(p*Urpp-0.75*(qpp0*UmB+2*qp*UmB_p+q*k*p*d2))/(0.75*(qppd4*UmB+q))
    return [d1,d2,d3,d4]

# test even: phi(0)=1,phi'(0)=0,phi'''(0)=0, vary phi''(0)=t; find A(pi)=0; check C(pi)
def AC_even(t):
    s=solve_ivp(rhs,[0,np.pi],[1,0,t,0],rtol=1e-11,atol=1e-13,method='Radau',max_step=0.005)
    phi,d1,d2,d3=s.y[:,-1]
    return k*phi+d2,-(k-1)*d1
print("even attempt (teammate III): t -> A(pi), C(pi)")
ts=np.linspace(-2,2,21);vals=[AC_even(t)[0] for t in ts]
for t in ts: print("  t=%.2f A=%.4f C=%.4f"%(t,AC_even(t)[0],AC_even(t)[1]))