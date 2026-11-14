import numpy as np, sys
sys.path.insert(0, '/work')
from crack3d import Operator2D, material, smooth_chi

mat = material()
Ep, mu = mat['Ep'], mat['mu']
s2p = np.sqrt(2*np.pi)
K10 = 1.0
C0 = 8*K10/(Ep*s2p)

Nx, Nz = 8192, 8
Lx = 8.0
dx = 2*Lx/Nx
xg = (np.arange(Nx)-Nx//2)*dx
zg = np.arange(Nz)*(1.0/Nz)
X, Z = np.meshgrid(xg, zg, indexing='ij')
op = Operator2D((Nx,Nz), dx, 1.0/Nz, Ep, mu)

# windowed straight-crack opening density (compact)
chi = smooth_chi(xg, -4.0, -6.0)     # 1 for x>-4, 0 for x<-6
w = np.where(xg<0, np.sqrt(np.maximum(-xg,0))*chi, 0.0)
by = (C0*w)[:,None]*np.ones((1,Nz))
bx = 0*by; bz = 0*by
tx,ty,tz = op.apply(bx,by,bz)
ty1 = ty[:,0].real

tinf = 1/np.sqrt(2*np.pi*np.abs(xg))
print("   x        t_y(grid)     -t_inf(x)     ratio   [x<0: expect ~1]")
for xq in [-3.0,-2.0,-1.0,-0.5,-0.2,-0.05,-0.01]:
    i = np.argmin(np.abs(xg-xq))
    print(f"{xq:8.2f} {ty1[i]:12.5f} {-tinf[i]:12.5f} {ty1[i]/(-tinf[i]):8.4f}")
print("   x>0: expect ~0")
for xq in [0.01,0.05,0.2,0.5,1.0,2.0]:
    i = np.argmin(np.abs(xg-xq))
    print(f"{xq:8.2f} {ty1[i]:12.5f} {0:12.5f}")
