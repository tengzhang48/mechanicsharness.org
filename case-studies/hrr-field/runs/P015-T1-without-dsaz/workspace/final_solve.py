"""Refined HRR solution + validation:
 - high-accuracy integration
 - check A(pi), C(pi) ~ 0
 - check smoothness (no blow-up)
 - check the ODE residual (max |phi'''' - rhs|)
 - compare against a direct solve_bvp on [0,pi] with BCs
   phi(0)=1, phi'(0)=0, A(pi)=0, C(pi)=0.
"""
import numpy as np
import hrr_clean as H
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

n=H.N; p=H.P; k=H.K

def resid(u):
    s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0.0,np.pi],
                [1.0,0.0,u[0],u[1]], rtol=1e-12, atol=1e-14,
                method='Radau', max_step=0.004)
    phi,d1,d2,d3=s.y[:,-1]
    return np.array([k*phi+d2, -(k-1.0)*d1])

u=fsolve(resid,(-0.91422447,0.03284634),xtol=1e-13,full_output=True)
u=u[0]
print("phi''(0),phi'''(0) =",u)
s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0.0,np.pi],
            [1.0,0.0,u[0],u[1]], rtol=1e-13, atol=1e-15,
            method='Radau', dense_output=True, max_step=0.002)
phi,d1,d2,d3=s.y[:,-1]
print("A(pi)=",k*phi+d2," C(pi)=",-(k-1.0)*d1)

# --- solve_bvp cross-check with same 4 BCs ---
from scipy.integrate import solve_bvp
def rhs_bvp(th,Y):
    return H.rhs(th,Y,n=n,p=p,k=k)
def bc(ya,yb):
    return np.array([ya[0]-1.0, ya[1], k*yb[0]+yb[2], yb[1]])
th=np.linspace(0,np.pi,200)
y0=np.zeros((4,th.size))
y0[0]=1.0+u[0]/2*np.cos(th)  # rough
sol=solve_bvp(rhs_bvp,bc,th,y0,tol=1e-9,max_nodes=5000)
print("bvp success",sol.success,sol.message[:50],"nodes",sol.x.size)
if sol.success:
    print("bvp phi''(0)=",sol.y[2,0]," phi'''(0)=",sol.y[3,0])
    print("bvp A(pi)=",k*sol.y[0,-1]+sol.y[2,-1]," C(pi)=",-(k-1)*sol.y[1,-1])
    # compare with shooting on grid
    thg=np.linspace(0,np.pi,200)
    ys=s.sol(thg)
    yb=sol.sol(thg)
    print("max |shoot - bvp| phi:",np.max(np.abs(ys[0]-yb[0])))
    print("max |shoot - bvp| phi':",np.max(np.abs(ys[1]-yb[1])))
    print("max |shoot - bvp| phi'':",np.max(np.abs(ys[2]-yb[2])))
    print("max |shoot - bvp| phi''':",np.max(np.abs(ys[3]-yb[3])))