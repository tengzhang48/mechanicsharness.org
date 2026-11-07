import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import leggauss
from explore import apply_R, Stilde

xg, wg = leggauss(2000); xg=0.5*(xg+1); wg=0.5*wg

def brute(P, Z):
    # integrate S(Z,xi) P(xi) dxi with singularity subtraction for the 1/(Z-xi) part
    # S = 1/(Z-xi) + Rreg ;  do PV of P/(Z-xi) by subtracting smooth part
    # PV int P/(Z-xi) = int (P(xi)-P(Z))/(Z-xi) dxi + P(Z)*PV int dxi/(Z-xi)
    val_s = np.sum(wg*(P(xg)-P(Z))/(Z-xg))
    val_s += P(Z)*np.log(Z/(1-Z))
    Rreg = lambda Z,X: (X**2-4*Z*X-Z**2)/(Z+X)**3
    val_r = np.sum(wg*Rreg(Z,xg)*P(xg))
    return -(val_s+val_r)

for P in [Poly([1.0]), Poly([0,1.0]), Poly([0,0,1.0]), Poly([1,0,1.0]), Poly([-1,3,-3,1.0])]:
    for Z in [0.1,0.3,0.5,0.7,0.9]:
        a=apply_R(P,Z,xg,wg); b=brute(P,Z)
        print("P=",np.round(P.coef,2),"Z=",Z," apply_R=%.8f brute=%.8f diff=%.1e"%(a,b,a-b))