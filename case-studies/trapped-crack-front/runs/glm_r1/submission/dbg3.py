import numpy as np, sys
sys.path.insert(0, '/work')
from crack3d import material, smooth_chi

mat = material(); Ep = mat['Ep']
Lcr = 2.0
C0y = 8*1.0/(Ep*np.sqrt(2*np.pi))

N = 4000
k = np.arange(N)
u = np.cos(np.pi*(k+0.5)/N)
sig = (u+1)*np.sqrt(Lcr)/2
x = -sig**2
chi = smooth_chi(x, -0.75*Lcr, -0.95*Lcr)
G = sig*chi*C0y

# direct Chebyshev projection
M = 64
c = np.zeros(M)
for m in range(M):
    Tm = np.cos(m*np.arccos(u))
    c[m] = (2.0/N)*np.sum(G*Tm)*(1.0 if m>0 else 0.5)
Gfit = np.polynomial.chebyshev.chebval(u, c)
print("proj err:", np.abs(Gfit-G).max())
m = np.arange(M)
Cex = (2/np.sqrt(Lcr))*np.sum(c*((-1.0)**(m-1))*m**2)
print(f"C extracted {Cex:.6f} vs C0 {C0y:.6f}")
print("coeffs:", np.round(c[:20],8))
# smoothness check: derivative of G wrt u near u=-1
print("G near u=-1:", G[:5]/C0y)
