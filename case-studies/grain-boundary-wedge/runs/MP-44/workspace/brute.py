import numpy as np
from numpy.polynomial.legendre import leggauss
import mpmath as mp
mp.mp.dps=30

def gA(k,xi):  # family A: (xi(1-xi))^{-1/2} T_k(2xi-1)
    xi=mp.mpf(xi)
    return (xi*(1-xi))**(-mp.mpf(1)/2)*mp.chebyt(k,2*xi-1)

def R_brute(k,Z):
    Z=mp.mpf(Z)
    # R(Z) = - [ PV int g/(Z-Xi) + int Rreg g ]
    # PV part via subtraction
    N=600
    x,w=leggauss(N)
    xi=0.5*(x+1.0); wt=0.5*w
    gv=np.array([float(gA(k,mp.mpf(v))) for v in xi])
    gZ=float(gA(k,mp.mpf(Z)))
    pv=np.sum(wt*(gv-gZ)/(float(Z)-xi))+gZ*np.log(abs(float(Z)/(1-float(Z))))
    Rgv=(xi**2-4*float(Z)*xi-float(Z)**2)/(float(Z)+xi)**3
    reg=np.sum(wt*Rgv*gv)
    return -(pv+reg)

for k in [1,2,3]:
    for Z in [0.2,0.5,0.7]:
        print("k=%d Z=%.1f  R_brute=%.6f"%(k,Z,R_brute(k,Z)))