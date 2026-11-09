"""2-parameter shooting: phi(0)=1, phi'(0)=0, free (phi''(0), phi'''(0)),
impose A(pi)=0 and C(pi)=0. Use fsolve on (phi''(0),phi'''(0)).
Try many initial guesses."""
import numpy as np, sys
import hrr_clean as H
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

n=H.N; p=H.P; k=H.K
def resid(u):
    s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0.0,np.pi],
                [1.0,0.0,u[0],u[1]], rtol=1e-10, atol=1e-13,
                method='Radau', max_step=0.01)
    phi,d1,d2,d3=s.y[:,-1]
    return np.array([k*phi+d2, -(k-1.0)*d1])

print("Solving 2-param system for p=1/(n+1)=",p)
sols=[]
for g in [(-1,0),(-0.5,0),(0,0),(0.5,0),(1,0),(-1,-1),(0,1),(0,-1),(1,1),(-0.7,0.3),(-0.6,-0.3)]:
    try:
        r=fsolve(resid, g, full_output=True)
        u=r[0]; info=r[1]; ok=r[2]
        f=resid(u)
        if np.linalg.norm(f)<1e-6:
            sols.append(tuple(np.round(u,6)))
            print(f"  guess {g} -> u={u} |resid|={np.linalg.norm(f):.2e} nfev={info['nfev']} id={id}")
    except Exception as e:
        print("  guess",g,"err",e)
    sys.stdout.flush()
print("unique solutions:",set(sols))