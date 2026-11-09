import numpy as np
from scipy.integrate import solve_ivp, quad
n=3.5; m=1/n; beta=1/(n+1); D=1+(beta-1)/n
def rhs(t,y):
    f,u=y; Q2=beta*beta*f*f+u*u; den=beta*beta*f*f+m*u*u
    return [u, -beta*f*((m-1)*beta*u*u+D*(Q2))/den]
sol=solve_ivp(rhs,[0,np.pi],[0,1],rtol=1e-13,atol=1e-15,method='DOP853',dense_output=True)
# closed form: rho(chi)=C0*(cos^2+m^2 sin^2)^(-1/(2(1+m))), Q=beta*rho
C0=(m**(1/(1+m)))/beta
def theta_of_chi(chi):
    g=lambda c:(np.cos(c)**2+m*np.sin(c)**2)/(np.cos(c)**2+m*m*np.sin(c)**2)
    val,_=quad(g,chi,np.pi/2)  # theta = (m+1)*int_chi^{pi/2} g dchi
    return (m+1)*val
def Q_closed(th):
    # invert theta(chi)=th by bisection on chi in [0,pi/2]
    lo,hi=0.0,np.pi/2
    for _ in range(80):
        mid=0.5*(lo+hi)
        if theta_of_chi(mid)>th: lo=mid
        else: hi=mid
    chi=0.5*(lo+hi)
    rho=C0*(np.cos(chi)**2+m*m*np.sin(chi)**2)**(-1/(2*(1+m)))
    return beta*rho
for deg in [0,30,60,90,120,150,180]:
    th=np.deg2rad(deg)
    Qn=sol.sol(th)[0]*0+np.sqrt(beta*beta*sol.sol(th)[0]**2+sol.sol(th)[1]**2)
    print(deg, Qn, Q_closed(th))
