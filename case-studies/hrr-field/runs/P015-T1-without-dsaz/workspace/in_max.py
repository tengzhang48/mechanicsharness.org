"""Try max-se=1 normalization and both signs; also compute Wbar."""
import numpy as np, hrr_corr as H
from scipy.integrate import solve_ivp
n=H.N;p=H.P;k=H.K
u=np.array([-0.9142244740,0.0328463401])
s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0,np.pi],[1,0,u[0],u[1]],
            rtol=1e-13,atol=1e-15,method='Radau',dense_output=True,max_step=0.002)
ths=np.linspace(1e-6,np.pi-1e-6,3000)
Ys=s.sol(ths); phi,d1,d2,d3=Ys
A=k*phi+d2;B=k*(k-1)*phi;C=-(k-1)*d1
se=np.sqrt(0.75*(A-B)**2+3*C**2)
semax=se.max()
print("se(0)=%.5f se_max=%.5f se(pi)=%.5f"%(se[0],semax,se[-1]))
def allq(th):
    Y=s.sol(th); phi,d1,d2,d3=Y
    A=k*phi+d2;B=k*(k-1)*phi;C=-(k-1)*d1
    se=np.sqrt(0.75*(A-B)**2+3*C**2)
    Ap=k*d1+d3;Bp=k*(k-1)*d1
    sep=(0.75*(A-B)*(Ap-Bp)+3*C*(-(k-1)*d2))/se
    q=se**(n-1);qp=(n-1)*se**(n-2)*sep
    Ur=0.75*q*(A-B)/p
    Urp=0.75*(qp*(A-B)+q*(Ap-Bp))/p
    Ut=(3*q*C-Urp)/(p-1)
    d4=H.rhs(th,Y)[3];App=k*d2+d4;Bpp=k*(k-1)*d2
    qpp=(n-1)*(n-2)*se**(n-3)*sep**2+(n-1)*se**(n-2)*((0.75*(2*(A-B)*(App-Bpp)+2*(Ap-Bp)**2)+3*(2*(-(k-1)*d2)**2+2*C*(-(k-1)*d3))-2*sep**2)/(2*se))
    Urpp=0.75*(qpp*(A-B)+2*qp*(Ap-Bp)+q*(App-Bpp))/p
    Utp=(3*(qp*C+q*C*(-(k-1)*d2))-Urpp)/(p-1) if False else (3*(qp*C+q*(-(k-1)*d2))-Urpp)/(p-1)
    return A,B,C,se,Ur,Urp,Urpp,Ut,Utp
# max-se=1
Sc=1.0/semax
def integ(trsign,Wfac):
    vals=[]
    for th in ths:
        A,B,C,se,Ur,Urp,Urpp,Ut,Utp=allq(th)
        a,b,c,sen=Sc*A,Sc*B,Sc*C,Sc*se
        urn,utn=Sc*Ur,Sc*Ut;urnp,utnp=Sc*Urp,Sc*Utp
        co,si=np.cos(th),np.sin(th)
        dur_dx=p*urn*co-urnp*si; dut_dx=p*utn*co-utnp*si
        dux_dx=dur_dx*co-dut_dx*si-(urn*si*si+utn*si*co)
        duy_dx=dur_dx*si+dut_dx*co+(urn*si*co-utn*si*si)
        tr=a*co+c*si; tt=c*co+b*si
        vals.append(Wfac*sen**(n+1.0)*co - trsign*(tr*dux_dx+tt*duy_dx))
    vals=np.array(vals)
    return 2*np.trapezoid(vals,ths)
print("I_n(+traction,n/(n+1))=",integ(+1,n/(n+1)))
print("I_n(-traction,n/(n+1))=",integ(-1,n/(n+1)))
# Wbar = int_{-pi}^{pi} Wtilde dth = 2 int_0^pi (n/(n+1)) se^{n+1} dth
Wbar=2*np.trapezoid((n/(n+1.0))*(Sc*se)**(n+1.0),ths)
print("Wbar=",Wbar)
print("B(0)/semax =",Sc*(k*(k-1)*1.0))