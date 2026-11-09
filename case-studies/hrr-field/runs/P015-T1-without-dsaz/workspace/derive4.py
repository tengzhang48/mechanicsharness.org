import numpy as np
a,b=1.3,0.7
def phi(t): return a*np.cos(1.5*t)+b*np.cos(0.5*t)
def d(tt,m):
    if m==0: return phi(tt)
    if m==1: return -1.5*a*np.sin(1.5*tt)-0.5*b*np.sin(0.5*tt)
    if m==2: return -2.25*a*np.cos(1.5*tt)-0.25*b*np.cos(0.5*tt)
    if m==3: return 3.375*a*np.sin(1.5*tt)+0.125*b*np.sin(0.5*tt)
p=0.5;k=1.5
for tt in [0.3,1.0,2.0,2.9]:
    pp=d(tt,0);d1=d(tt,1);d2=d(tt,2);d3=d(tt,3)
    A=k*pp+d2;B=k*(k-1)*pp;C=-(k-1)*d1
    Ap=k*d1+d3;Bp=k*(k-1)*d1;Cp=-(k-1)*d2
    UmB=A-B;UmB_p=Ap-Bp
    Ur=UmB/(2*p);Urp=UmB_p/(2*p)
    Ut=(2*C-Urp)/(p-1)
    Utp_II=(B-A)/2-Ur
    h=1e-5
    pp2=d(tt+h,0);d12=d(tt+h,1);d22=d(tt+h,2);d32=d(tt+h,3)
    A2=k*pp2+d22;B2=k*(k-1)*pp2;UmB2=A2-B2
    Ap2=k*d12+d32;Bp2=k*(k-1)*d12
    Urp2=(Ap2-Bp2)/(2*p)
    Urpp=(Urp2-Urp)/h
    Cpp=-(k-1)*d22
    Utp_diff=(2*Cp-Urpp)/(p-1)
    print(f"t={tt}: IIt={Utp_II:.6f} difft={Utp_diff:.6f} diff={Utp_II-Utp_diff:.2e}")
# Also just integrate Ut numerically and compare
from scipy.integrate import cumulative_trapezoid
ts=np.linspace(0,3.2,20001)
def UT(tt):
    pp=d(tt,0);d1=d(tt,1);d2=d(tt,2);d3=d(tt,3)
    A=k*pp+d2;B=k*(k-1)*pp;C=-(k-1)*d1
    Ap=k*d1+d3;Bp=k*(k-1)*d1
    Ur=(A-B)/(2*p);Urp=(Ap-Bp)/(2*p)
    return Ur,(B-A)/2-Ur
vals=np.array([UT(t) for t in ts])
# Ut' target from II = vals[:,1]; integrate -> Ut_II
UtII=cumulative_trapezoid(vals[:,1],ts,initial=0)
# Ut from III = (2C-Urp)/(p-1)
UtIII=np.array([ (2*(-(k-1)*d(t,1)) - (k* d(t,1)+d(t,3)-(k*(k-1))*d(t,1))/(2*p))/(p-1) for t in ts])
print("max |Ut_II - Ut_III| over range:",np.max(np.abs(UtII-UtIII)))