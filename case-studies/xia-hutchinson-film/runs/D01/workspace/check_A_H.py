import numpy as np
from j2d import solve, Jarc
nu=0.3; Ebar=1.0; s0=1.0; ell=1.0; h=1.0
Lx=4.0; Nx=600
print("CASE A H-dependence (ell=1): analytic gamma_A=tanh(H/(2ell))")
for H in [0.5,1.0,2.0]:
    Wd=H; Ny=100; dy=Wd/Ny; x2=np.arange(Ny+1)*dy
    Bw=-s0*ell/(Ebar*np.cosh(H/(2*ell)))
    left=Bw*np.sinh((x2-H/2)/ell)
    U,dx,dyy=solve(Wd,Lx,Nx,Ny,left,np.zeros(Ny+1),'mix','mix')
    GA=(s0**2*ell/Ebar)*np.tanh(H/(2*ell))
    r0=0.25*min(H,ell*2)
    Gnum=2*Jarc(U,dx,dyy,r0,H/2 if False else 0.0)/h  # tip at x2=0
    print(f"H={H}: r0={r0:.3f} G_num={Gnum:.6e} analytic={GA:.6e} rel.err={abs(Gnum-GA)/GA:.2e}")
