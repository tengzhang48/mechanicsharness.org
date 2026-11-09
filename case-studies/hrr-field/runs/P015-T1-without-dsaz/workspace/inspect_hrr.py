"""Build the HRR Mode I solution for p=1/(n+1), inspect symmetry & shape."""
import numpy as np, sys
import hrr_clean as H
from scipy.integrate import solve_ivp

n=H.N; p=H.P; k=H.K
u0=(-0.91422447,0.03284634)
s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0.0,np.pi],
            [1.0,0.0,u0[0],u0[1]], rtol=1e-12, atol=1e-14,
            method='Radau', dense_output=True, max_step=0.005)
th=np.linspace(0,np.pi,19)
Y=s.sol(th)
phi=Y[0]; d1=Y[1]; d2=Y[2]; d3=Y[3]
A=k*phi+d2; B=k*(k-1)*phi; C=-(k-1)*d1
se=np.sqrt(0.75*(A-B)**2+3*C**2)
print(" theta     phi        A          B          C          se")
for i in range(len(th)):
    print(f" {th[i]:.3f}  {phi[i]:+.5f}  {A[i]:+.5f}  {B[i]:+.5f}  {C[i]:+.5f}  {se[i]:+.5f}")
print("A(pi)=",k*phi[-1]+d2[-1]," C(pi)=",-(k-1)*d1[-1])
print("phi'''(0)=",u0[1]," d1(0)=",d1[0]," C(0)=",-(k-1)*d1[0])
# check evenness: extend to negative theta using ODE from 0 with -theta?
# Instead check phi'(theta) values: even sol would have phi'(pi)=? not informative.
# Check sigma_rt at 0.
print("se(0)=",se[0]," se at pi=",se[-1])