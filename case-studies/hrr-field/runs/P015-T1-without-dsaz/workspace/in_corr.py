"""Corrected I_n: U_r=(3/4)q(A-B)/p, U_t=(3qC-U_r')/(p-1).
Also try the overall-const scaling. Check sign of I_n."""
import numpy as np, hrr_corr as H
from scipy.integrate import solve_ivp
n=H.N;p=H.P;k=H.K
u=np.array([-0.9142244740,0.0328463401])
s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0,np.pi],[1,0,u[0],u[1]],
            rtol=1e-13,atol=1e-15,method='Radau',dense_output=True,max_step=0.002)
def allq(th):
    Y=s.sol(th); phi,d1,d2,d3=Y
    A=k*phi+d2;B=k*(k-1)*phi;C=-(k-1)*d1
    se=np.sqrt(0.75*(A-B)**2+3*C**2)
    Ap=k*d1+d3;Bp=k*(k-1)*d1;Cp=-(k-1)*d2
    sep=(0.75*(A-B)*(Ap-Bp)+3*C*Cp)/se
    d4=H.rhs(th,Y)[3];App=k*d2+d4;Bpp=k*(k-1)*d2
    q=se**(n-1);qp=(n-1)*se**(n-2)*sep
    qpp=(n-1)*(n-2)*se**(n-3)*sep**2+(n-1)*se**(n-2)*((0.75*(2*(A-B)*(App-Bpp)+2*(Ap-Bp)**2)+3*(2*Cp*Cp+2*C*(-(k-1)*d3))-2*sep**2)/(2*se))
    Ur=0.75*q*(A-B)/p
    Urp=0.75*(qp*(A-B)+q*(Ap-Bp))/p
    Urpp=0.75*(qpp*(A-B)+2*qp*(Ap-Bp)+q*(App-Bpp))/p
    Ut=(3*q*C-Urp)/(p-1)
    Utp=(3*(qp*C+q*Cp)-Urpp)/(p-1)
    return A,B,C,se,Ur,Urp,Urpp,Ut,Utp
A0,B0,C0,se0,_,_,_,_,_=allq(0.0);Sc=1.0/se0
ths=np.linspace(1e-6,np.pi-1e-6,3000)
Wl=[];Tl=[]
for th in ths:
    A,B,C,se,Ur,Urp,Urpp,Ut,Utp=allq(th)
    a,b,c,sen=Sc*A,Sc*B,Sc*C,Sc*se
    urn,utn=Sc*Ur,Sc*Ut;urnp,utnp=Sc*Urp,Sc*Utp
    co,si=np.cos(th),np.sin(th)
    dur_dx=p*urn*co-urnp*si; dut_dx=p*utn*co-utnp*si
    dux_dx=dur_dx*co-dut_dx*si-(urn*si*si+utn*si*co)
    duy_dx=dur_dx*si+dut_dx*co+(urn*si*co-utn*si*si)
    tr=a*co+c*si; tt=c*co+b*si
    Wl.append((n/(n+1.0))*sen**(n+1.0)*co)
    Tl.append(tr*dux_dx+tt*duy_dx)
Wl=np.array(Wl);Tl=np.array(Tl)
iW=np.trapezoid(Wl,ths); iT=np.trapezoid(Tl,ths)
print("intW[0,pi]=%.5f intT[0,pi]=%.5f"%(iW,iT))
print("I_n full=%.5f"%(2*(iW-iT)))
print("A(0)/se0=%.5f B(0)/se0=%.5f se0=%.6f"%(A0/se0,B0/se0,se0))