"""Griffith crack (half-length L=1, unit pressure): check
   int sigma_kk dA = -pi L^2 (large disk)   and   int sigma_kk^2 dA = 2 pi L^2,
used in the reciprocal-theorem check of V2(inf) = -pi*eta*(1-nu)."""
import numpy as np
from scipy import integrate
def integrand(xi, ph, power):
    w = xi + 1j*ph
    Z = 1/np.tanh(w) - 1.0
    skk = 2*np.real(Z)
    jac = np.abs(np.sinh(w))**2
    return skk**power*jac
I2, e2 = integrate.dblquad(lambda ph, xi: integrand(xi, ph, 2), 0, 40, 0, 2*np.pi, epsabs=1e-11, epsrel=1e-11)
print('int skk^2 dA =', I2, ' 2*pi =', 2*np.pi)
for Xi in [5.0, 10.0]:
    I1, e1 = integrate.dblquad(lambda ph, xi: integrand(xi, ph, 1), 0, Xi, 0, 2*np.pi, epsabs=1e-11, epsrel=1e-11)
    print('int skk dA over ellipse xi<%g =' % Xi, I1, ' -pi =', -np.pi)
