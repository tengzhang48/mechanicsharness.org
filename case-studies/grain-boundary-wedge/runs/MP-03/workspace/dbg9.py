import numpy as np
from eig_g import tquad
from gbcore import Rfun
from numpy.polynomial.legendre import leggauss

Nv, Mt = 8, 400
xc = 0.5*(1-np.cos(np.pi*np.arange(1,Nv+2)/(Nv+1)))
t, wt = tquad(0.0,1.0,Mt); xi = 1-t**2
x = xc[5]
den = x-1+t**2
print("min |den|:", np.min(np.abs(den)), " num zeros near?", )
# check exact-zero den
print("any den==0:", np.any(den==0))
