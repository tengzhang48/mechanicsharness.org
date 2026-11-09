import sys, numpy as np
sys.path.insert(0, 'knowledge/checks')
from crack_tip_checks import check_field, format_report

E=70000.0; nu=0.30
# effective nu from solution (T solved below recomputed here quickly)
from scipy.optimize import brentq
Czz=-0.005; J=0.8
def nu_eff_of(T):
    Cdev=(1+nu)/E+1.5*(1/T-1/E); Cbulk=(1-2*nu)/E
    G=1/(2*Cdev); Kb=1/(3*Cbulk)
    return (3*Kb-2*G)/(2*(3*Kb+G))
def F(T):
    Eeff=1/((2/3)*((1+nu)/E+1.5*(1/T-1/E))+(1/3)*((1-2*nu)/E))
    Atr=nu/E+0.5*(1/T-1/E)
    return Atr**2*Eeff - Czz**2*np.pi/(2*J)
T=brentq(F,1,E-1,xtol=1e-14)
nu_eff=nu_eff_of(T)
print("T",T,"nu_eff",nu_eff)

th=np.linspace(0,np.pi,20001)
c=np.cos(th/2); s=np.sin(th/2); c3=np.cos(3*th/2); s3=np.sin(3*th/2)
Srr=5/4*c-1/4*c3
Stt=3/4*c+1/4*c3
Srt=1/4*s+1/4*s3
rep=check_field(th,Srr,Stt,Srt,p=0.5,n=1.0,beta=np.pi,plane="stress",
                law="linear",nu=nu_eff,mode="I",leading=True,
                reported_J=None)
print(format_report(rep))
