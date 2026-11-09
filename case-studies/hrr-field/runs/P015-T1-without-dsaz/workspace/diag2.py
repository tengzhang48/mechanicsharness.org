import numpy as np, hrr_clean as H
from scipy.integrate import solve_ivp
n=H.N;p=H.P;k=H.K
u=np.array([-0.9142244740,0.0328463401])
s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0,np.pi],[1,0,u[0],u[1]],
            rtol=1e-13,atol=1e-15,method='Radau',dense_output=True,max_step=0.002)
def UrUt(th,Y=None):
    phi,d1,d2,d3=s.sol(th) if Y is None else Y
    A=k*phi+d2; B=k*(k-1)*phi; C=-(k-1)*d1
    se=np.sqrt(0.75*(A-B)**2+3*C**2)
    return A,B,C,se,d1,d2,d3
# Two ways to get Ut':
def Ut_from_III(th):
    A,B,C,se,d1,d2,d3=UrUt(th)
    q=se**(n-1)
    Ap=k*d1+d3;Bp=k*(k-1)*d1;Cp=-(k-1)*d2
    sep=(0.75*(A-B)*(Ap-Bp)+3*C*Cp)/se
    qp=(n-1)*se**(n-2)*sep
    Urp=(qp*(A-B)+q*(Ap-Bp))/(2*p)
    Ut=(2*q*C-Urp)/(p-1)
    return Ut
def Ut_from_III_prime(th):
    # d/dth of Ut above; compute Ur'' from RHS
    h=1e-6
    return (Ut_from_III(th+h)-Ut_from_III(th-h))/(2*h)
def Utp_from_II(th):
    A,B,C,se,d1,d2,d3=UrUt(th)
    q=se**(n-1)
    Ur=q*(A-B)/(2*p)
    return q*(B-A)/2 - Ur
for th in [0.5,1.2,2.0,2.8]:
    print(f"th={th}: Ut'_III={Ut_from_III_prime(th):.5f}  Ut'_II={Utp_from_II(th):.5f}  diff={Ut_from_III_prime(th)-Utp_from_II(th):.2e}")
# also verify ODE by finite differencing phi''' from rhs
print("ODE self-consistency: compare d3' (rhs) vs fd d3:")
h=1e-6
for th in [0.5,1.2,2.0,2.8]:
    Y=s.sol(th); Yp=s.sol(th+h); Ym=s.sol(th-h)
    d3p_fd=(Yp[3]-Ym[3])/(2*h)
    d3p_rhs=H.rhs(th,Y)[3]
    print(f"   th={th}: d3'_fd={d3p_fd:.5f} d3'_rhs={d3p_rhs:.5f} diff={d3p_fd-d3p_rhs:.2e}")