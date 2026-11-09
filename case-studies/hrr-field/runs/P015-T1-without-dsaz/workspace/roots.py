"""Find ALL roots of the 2x2 system (phi''(0),phi'''(0)) via many fsolve seeds."""
import numpy as np, hrr_clean as H
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
n=H.N;p=H.P;k=H.K
def resid(u):
    s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0.0,np.pi],
                [1.0,0.0,u[0],u[1]], rtol=1e-9, atol=1e-11,
                method='Radau', max_step=0.01)
    if not s.success: return np.array([1e3,1e3])
    phi,d1,d2,d3=s.y[:,-1]
    return np.array([k*phi+d2, -(k-1.0)*d1])
roots=[]
for a in np.linspace(-3,1.5,10):
    for b in np.linspace(-1,1,5):
        try:
            sol=fsolve(resid,(a,b),xtol=1e-10,full_output=True)
            u=sol[0]; ier=sol[2]
            r=resid(u)
            if np.linalg.norm(r)<1e-6:
                # dedupe
                if not any(np.allclose(u,z,atol=1e-4) for z in roots):
                    roots.append(u)
        except Exception: pass
print("distinct roots (phi''(0),phi'''(0)):")
for z in roots:
    # evaluate B(0),se(0) and shape
    s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0,np.pi],
                [1.0,0.0,z[0],z[1]], rtol=1e-11, atol=1e-13, method='Radau', max_step=0.005)
    phi,d1,d2,d3=s.sol(0.0)
    A=k*phi+d2;B=k*(k-1.0)*phi;C=-(k-1.0)*d1
    se=np.sqrt(0.75*(A-B)**2+3*C**2)
    # check se>0 everywhere and C(0)=0
    ths=np.linspace(0,np.pi,100); Ys=s.sol(ths)
    Aa=k*Ys[0]+Ys[2];Bb=k*(k-1)*Ys[0];Cc=-(k-1)*Ys[1]
    ses=np.sqrt(0.75*(Aa-Bb)**2+3*Cc**2)
    print(f"  u=({z[0]:+.5f},{z[1]:+.5f}) B(0)/se(0)={B/se:.5f} A(0)/se(0)={A/se:.5f} se_min={ses.min():.3e} se_max={ses.max():.3f} A(pi)={Aa[-1]:.1e} C(pi)={Cc[-1]:.1e}")