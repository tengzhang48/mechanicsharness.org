import numpy as np, hrr_clean as H
from scipy.integrate import solve_ivp
n=H.N;p=H.P;k=H.K
u=np.array([-0.9142244740,0.0328463401])
s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0,np.pi],[1,0,u[0],u[1]],
            rtol=1e-13,atol=1e-15,method='Radau',dense_output=True,max_step=0.002)
def parts(th):
    phi,d1,d2,d3=s.sol(th)
    A=k*phi+d2;B=k*(k-1)*phi;C=-(k-1)*d1
    se=np.sqrt(0.75*(A-B)**2+3*C**2)
    Ap=k*d1+d3;Bp=k*(k-1)*d1;Cp=-(k-1)*d2
    d4=H.rhs(th,s.sol(th))[3]
    App=k*d2+d4;Bpp=k*(k-1)*d2;Cpp=-(k-1)*d3
    sep=(0.75*(A-B)*(Ap-Bp)+3*C*Cp)/se
    return A,B,C,se,Ap,Bp,Cp,App,Bpp,Cpp,sep
for th in [0.5,1.2,2.0,2.8]:
    A,B,C,se,Ap,Bp,Cp,App,Bpp,Cpp,sep=parts(th)
    q=se**(n-1); qp=(n-1)*se**(n-2)*sep
    # se''
    se2pp=0.75*(2*Ap*Ap+2*(A-B)*(App) - 2*Bp*Bp- ... ) if False else None
    qpp=None
    Ur=q*(A-B)/(2*p)
    Urp=(qp*(A-B)+q*(Ap-Bp))/(2*p)
    # Ur'' via ODE form
    Urpp=2*(qp*C+q*Cp)-(p-1)*(q*(B-A)/2-Ur)
    # Ut' from III differentiated: [2(qp C+ q Cp)-Urpp]/(p-1)
    Utp_III=(2*(qp*C+q*Cp)-Urpp)/(p-1)
    Utp_II=q*(B-A)/2-Ur
    print(f"th={th}: Ut'_III={Utp_III:.6f} Ut'_II={Utp_II:.6f} diff={Utp_III-Utp_II:.2e}")