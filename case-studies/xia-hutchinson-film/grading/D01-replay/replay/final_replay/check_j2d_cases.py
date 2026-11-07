import numpy as np
from j2d import solve, Jarc

nu=0.3; Ebar=1.0; s0=1.0; ell=1.0; h=1.0

# ---------------- Case A: verify G_A ----------------
H=1.0; Wd=H; Lx=4.0; Nx=400; Ny=50
x2A=np.arange(Ny+1)*(Wd/Ny)
Bw=-s0*ell/(Ebar*np.cosh(H/(2*ell)))
wakeA=Bw*np.sinh((x2A-H/2)/ell)
UA,dxA,dyA=solve(Wd,Lx,Nx,Ny,wakeA,np.zeros(Ny+1),'mix','mix')
GA=(s0**2*ell/Ebar)*np.tanh(H/(2*ell))
print("CASE A: analytic G_A =",GA)
for r0 in [0.2,0.3,0.4,0.5]:
    J=Jarc(UA,dxA,dyA,r0,0.0)
    Gnum=2*J/h
    print(f"  r0={r0:3.1f}  G_num=2*J_arc={Gnum:.6e}  rel.err={abs(Gnum-GA)/GA:.2e}")

# ---------------- Case B: verify G_B ----------------
# shifted domain y=x2-H in [0,H]; new crack at y=0, old crack at y=H
Wd=H; Ny=100; Nx=600; Lx=4.0
dy=Wd/Ny
x2B=np.arange(Ny+1)*dy
Ba=-s0*ell/(Ebar*np.cosh(H/ell))      # old-crack wake (strip width 2H), midplane y=0
Bb=-s0*ell/(Ebar*np.cosh(H/(2*ell)))  # split wake (strip width H)
right_u2=Ba*np.sinh(x2B/ell)          # far ahead: old-crack wake
left_u2 =Bb*np.sinh((x2B-H/2)/ell)    # far behind: split strips
UB,dxB,dyB=solve(Wd,Lx,Nx,Ny,left_u2,right_u2,'mix','free')
GB=(s0**2*ell/Ebar)*(2*np.tanh(H/(2*ell))-np.tanh(H/ell))
print("CASE B: analytic G_B =",GB)
for r0 in [0.15,0.2,0.25,0.3]:
    J=Jarc(UB,dxB,dyB,r0,0.0)
    Gnum=2*J/h
    print(f"  r0={r0:4.2f}  G_num=2*J_arc={Gnum:.6e}  rel.err={abs(Gnum-GB)/GB:.2e}")
