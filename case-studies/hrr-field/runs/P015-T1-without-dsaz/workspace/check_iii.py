"""Test which (III) is consistent: integrate U_t' from (II) and compare to the
algebraic U_t from (III) for BOTH sign conventions, using my phi-field."""
import numpy as np, hrr_clean as C
from scipy.integrate import solve_ivp, cumulative_trapezoid
n,p,k=C.N,C.P,C.K
u=np.array([-0.9142244740,0.0328463401])
s=solve_ivp(lambda th,Y: C.rhs(th,Y,n=n,p=p,k=k),[0,np.pi],[1,0,u[0],u[1]],
            rtol=1e-13,atol=1e-15,method='Radau',dense_output=True,max_step=0.001)
ths=np.linspace(0,np.pi,5000)
Y=s.sol(ths);phi,d1,d2,d3=Y
A=k*phi+d2;B=k*(k-1)*phi;Cc=-(k-1)*d1
se=np.sqrt(0.75*(A-B)**2+3*Cc**2)
Ap=k*d1+d3;Bp=k*(k-1)*d1;Cp=-(k-1)*d2
sep=(0.75*(A-B)*(Ap-Bp)+3*Cc*Cp)/se
q=se**(n-1);qp=(n-1)*se**(n-2)*sep
Ur=q*(A-B)/(2*p)          # using the 1/2 form (scale irrelevant for consistency)
Urp=(qp*(A-B)+q*(Ap-Bp))/(2*p)
Utp_II=q*(B-A)/2 - Ur
# candidate U_t from (III)-mine: Ur'+(p-1)Ut=2qC  -> Ut=(2qC-Ur')/(p-1)
Ut_mine=(2*q*Cc-Urp)/(p-1)
Utp_mine=(2*(qp*Cc+q*Cp)-( (n-1)*(n-2)*se**(n-3)*sep**2+(n-1)*se**(n-2)*(0)/1 ))/1 if False else None
# candidate U_t from (III)-teammate: Ur'-(p+1)Ut=2qC -> Ut=(Ur'-2qC)/(p+1)
Ut_team=(Urp-2*q*Cc)/(p+1)
# integrate Utp_II to get U_t^II up to a constant; compare SHAPE via derivative
Ut_II=cumulative_trapezoid(Utp_II,ths,initial=0.0)
# compare d/dth(Ut_mine) to Utp_II
dUt_mine=np.gradient(Ut_mine,ths)
dUt_team=np.gradient(Ut_team,ths)
print("checking (Ut)' == Utp_II:")
for i in [700,1500,2500,3500,4500]:
    print("  th=%.3f  Utp_II=%.5f  dUt_mine=%.5f  dUt_team=%.5f"%(ths[i],Utp_II[i],dUt_mine[i],dUt_team[i]))
print("max|dUt_mine-Utp_II| =",np.max(np.abs(dUt_mine-Utp_II)))
print("max|dUt_team-Utp_II| =",np.max(np.abs(dUt_team-Utp_II)))