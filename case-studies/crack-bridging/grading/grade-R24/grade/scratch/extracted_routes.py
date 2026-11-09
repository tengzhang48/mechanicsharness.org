import numpy as np

from scipy.optimize import brentq

from scipy.integrate import quad

E    = 130e9

nu   = 0.28

Ep   = E/(1-nu**2)

a    = 2.40e-3

t    = 160e6

dc   = 9.00e-6

G0   = 22.0

K0   = np.sqrt(Ep*G0)

def Sigma_of_c(c):           # K_tip = K0 condition
    return K0/np.sqrt(np.pi*c) + (2*t/np.pi)*np.arccos(a/c)

def delta_t_A(L):            # route A closed form
    c = a + L
    bracket = K0/np.sqrt(np.pi*c)                      # = Sigma - (2t/pi)arccos
    return (4*bracket/Ep)*np.sqrt(c**2 - a**2) + (8*t*a/(np.pi*Ep))*np.log(1 + L/a)

def delta_t_B(L):            # route B: ln-kernel quadrature, one adaptive call
    c = a + L
    S = Sigma_of_c(c)
    Rx = np.sqrt(c**2 - a**2)
    def G(s):
        Rs = np.sqrt(c**2 - s**2)
        return np.log((Rx + Rs)/abs(Rx - Rs))
    # log singularity at s=a (integrable), G->0 at s=c; split to help quad
    tot = quad(G, a, c, points=[a, 0.5*(a+c)], limit=100,
               epsabs=1e-13, epsrel=1e-12)[0]
    # remote part over whole slit integrates to 4*Sigma*sqrt(c^2-a^2)/Ep (exact)
    return (4*S/Ep)*Rx - (4*t/(np.pi*Ep))*tot

LA = brentq(lambda L: delta_t_A(L) - dc, 1e-6, 0.2, xtol=1e-14, rtol=1e-14)

LB = brentq(lambda L: delta_t_B(L) - dc, 1e-6, 0.2, xtol=1e-14, rtol=1e-14)
