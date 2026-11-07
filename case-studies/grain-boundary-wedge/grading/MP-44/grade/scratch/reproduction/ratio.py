import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import leggauss
from explore import apply_R

xg,wg=leggauss(2000); xg=0.5*(xg+1); wg=0.5*wg
base=Poly([-1.0,3,-3,1.0])  # (Z-1)^3
for p in [1,2,3,4,5]:
    f=base*Poly([0]*p+[1.0])
    Zt=np.linspace(0.05,0.95,13)
    Rf=np.array([apply_R(f.deriv(3),Z,xg,wg) for Z in Zt])
    fv=np.array([f(Z) for Z in Zt])
    print("f=Z^%d(Z-1)^3  R[f]/f:"%p, np.round(Rf/fv,3))