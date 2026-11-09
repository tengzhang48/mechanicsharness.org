"""Test I_n conventions: sign of traction term, and the (n/(n+1)) factor."""
import numpy as np, hrr_clean as H
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
    q=se**(n-1); qp=(n-1)*se**(n-2)*sep
    Ur=q*(A-B)/(2*p)
    Urp=(qp*(A-B)+q*(Ap-Bp))/(2*p)
    Ut=(2*q*C-Urp)/(p-1)
    Utp=(2*(qp*C+q*Cp)-(2*(qp*(A-B)+q*(Ap-Bp))/(2*p) if False else 0))/1
    # Utp via ODE:
    d4=H.rhs(th,Y)[3]; App=k*d2+d4;Bpp=k*(k-1)*d2
    Urpp=( (n-1)*(n-2)*se**(n-3)*sep**2+(n-1)*se**(n-2)*((0.75*(2*(A-B)*(App-Bpp)+2*(Ap-Bp)**2)+3*(2*Cp*Cp+2*C*(-(k-1)*d3))-2*sep**2)/(2*se)) )*(A-B)
    Urpp=(Urpp+2*qp*(Ap-Bp)+q*(App-Bpp))/(2*p)
    Utp=(2*(qp*C+q*Cp)-Urpp)/(p-1)
    return A,B,C,se,Ur,Urp,Ut,Utp
A0,B0,C0,se0,_,_,_,_=allq(0.0); Sc=1.0/se0
ths=np.linspace(1e-6,np.pi-1e-6,3000)
def build(sign_tr, Wfac):
    vals=[]
    for th in ths:
        A,B,C,se,Ur,Urp,Ut,Utp=allq(th)
        a,b,c,sen=Sc*A,Sc*B,Sc*C,Sc*se
        urn,utn=Sc*Ur,Sc*Ut; urnp,utnp=Sc*Urp,Sc*Utp
        dUr_dx=p*urn*np.cos(th)-urnp*np.sin(th)
        dUt_dx=p*utn*np.cos(th)-utnp*np.sin(th)
        dux_dx=dUr_dx*np.cos(th)-dUt_dx*np.sin(th)
        duy_dx=dUr_dx*np.sin(th)+dUt_dx*np.cos(th)
        tr=a*np.cos(th)+c*np.sin(th); tt=c*np.cos(th)+b*np.sin(th)
        Wt=Wfac*sen**(n+1.0)*np.cos(th)
        vals.append(Wt-sign_tr*(tr*dux_dx+tt*duy_dx))
    return 2*np.trapezoid(np.array(vals),ths)
print("I_n(+traction, n/(n+1)) =",build(+1,n/(n+1)))
print("I_n(-traction, n/(n+1)) =",build(-1,n/(n+1)))
print("I_n(+traction, 1)       =",build(+1,1.0))
print("I_n(-traction, 1)       =",build(-1,1.0))