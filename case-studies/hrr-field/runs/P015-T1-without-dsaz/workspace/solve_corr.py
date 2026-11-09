"""Solve corrected HRR (factor 3/2). 2 free ICs, 2 free-surface BCs."""
import numpy as np, hrr_corr as H
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
n=H.N;p=H.P;k=H.K
def resid(u):
    s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0.0,np.pi],
                [1.0,0.0,u[0],u[1]], rtol=1e-11, atol=1e-13,
                method='Radau', max_step=0.005)
    if not s.success: return np.array([1e3,1e3])
    phi,d1,d2,d3=s.y[:,-1]
    return np.array([k*phi+d2, -(k-1.0)*d1])
print("seeds -> roots:")
for seed in [(-0.7050162,0.0932568),(-0.91422447,0.03284634),(-0.7,0.1),(0.7,0.0),(-1.8,1.5)]:
    sol=fsolve(resid,seed,xtol=1e-13,full_output=True)
    u=sol[0]; r=resid(u)
    print(f"  seed {seed} -> u=({u[0]:+.6f},{u[1]:+.6f}) |r|={np.linalg.norm(r):.2e}")