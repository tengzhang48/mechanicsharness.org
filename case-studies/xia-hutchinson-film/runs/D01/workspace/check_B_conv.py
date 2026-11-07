import numpy as np
from j2d import solve, Jarc

nu=0.3; Ebar=1.0; s0=1.0; ell=1.0; h=1.0
H=1.0; Wd=H
GB=(s0**2*ell/Ebar)*(2*np.tanh(H/(2*ell))-np.tanh(H/ell))
print("analytic G_B =",GB)
for (Lx,Nx,Ny) in [(2.5,250,50),(2.5,500,100),(2.5,750,150)]:
    dy=Wd/Ny; x2B=np.arange(Ny+1)*dy
    Ba=-s0*ell/(Ebar*np.cosh(H/ell)); Bb=-s0*ell/(Ebar*np.cosh(H/(2*ell)))
    right=Ba*np.sinh(x2B/ell); left=Bb*np.sinh((x2B-H/2)/ell)
    UB,dx,dyy=solve(Wd,Lx,Nx,Ny,left,right,'mix','free')
    for r0 in [0.25,0.35]:
        Gnum=2*Jarc(UB,dx,dyy,r0,0.0)/h
        print(f"Lx={Lx} Nx={Nx} Ny={Ny} dx={dx:.4f} r0={r0}: G_num={Gnum:.6e} rel.err={abs(Gnum-GB)/GB:.3e}")
