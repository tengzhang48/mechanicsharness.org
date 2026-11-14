import numpy as np, sys
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem, material, smooth_chi

mat = material()
p = WavyFrontProblem(A=0.0, harm=1, Nxf=4096, Nz=64, M=16, Nzmode=4,
                     Lcr=2.0, Lx=3.5)
Ep = p.Ep; s2p = np.sqrt(2*np.pi)
C0y = 8*p.K10/(Ep*s2p)

# exact G at a=0: G(x) = C0*sqrt(-x)*chi(x) for x<0
chi = smooth_chi(p.xg, -0.75*p.Lcr, -0.95*p.Lcr)
Gex = np.where(p.xg < 0, np.sqrt(np.maximum(-p.xg,0))*chi, 0.0) * C0y

# project onto basis (component y, z-mode 0): values on grid
# u(x) = 2*sqrt(-x)/sqrt(Lcr)-1
sig = np.sqrt(np.maximum(-p.xg, 0.0))
u = np.clip(2*sig/np.sqrt(p.Lcr)-1, -1, 1)
sel = (p.xg < 0) & (p.xg > -p.Lcr)
# Chebyshev fit of Gex/C0 in u on sel
from numpy.polynomial.chebyshev import chebfit, chebval
ccoef = chebfit(u[sel], (Gex[sel]/1.0), p.M-1)
print("cheb coeffs (first 8):", np.round(ccoef[:8],5))

# extraction: C = (2/sqrt(Lcr)) * sum_m c_m * (-1)^{m-1} m^2
m = np.arange(p.M)
Cex = (2/np.sqrt(p.Lcr)) * np.sum(ccoef*((-1.0)**(m-1))*m**2)
print(f"extracted C = {Cex:.6f}   exact C0_y = {C0y:.6f}  ratio {Cex/C0y:.4f}")

# also check Chebyshev fit quality at u=-1 region
fitv = chebval(u, ccoef)
err = np.abs(fitv[sel]-Gex[sel]).max()
print(f"fit max err {err:.2e}, Gex max {Gex.max():.3f}")
