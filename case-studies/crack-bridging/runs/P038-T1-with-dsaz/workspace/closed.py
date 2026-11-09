"""Exact closed-form solve using the partially-bridged crack solution
(eqs 45,46 of Bao & Suo 1992; Tada et al.).

Uniform closure sigma0 on a<|x|<eta*a (eta = 1 + L/a), remote sigma_bar.
  delta_t E'/(s0 a) = 4 sqrt(eta^2-1)[ s/s0 - (2/pi) arccos(1/eta) ] + (8/pi) ln eta   (45)
  K_tip/(s0 sqrt(pi a)) = sqrt(eta)[ s/s0 - (2/pi) arccos(1/eta) ]                        (46)

Conditions at first connector failure (root |x|=a):
  delta_t = delta_c        (connector adjacent to slit reaches delta_c)
  K_tip   = K0 = sqrt(E' Gamma0)
"""
import numpy as np
from scipy.optimize import brentq

E=130e9; nu=0.28; Ep=E/(1-nu**2); s0=160e6; dc=9e-6; G0=22.0; a=2.4e-3
K0=np.sqrt(Ep*G0)
ell=Ep*dc/s0

def bracket(eta):
    r = arccos_term = (2/np.pi)*np.arccos(1/eta)
    return r

def sbar_from_ktip(eta):
    # (46): K_tip = s0 sqrt(pi a) sqrt(eta)[ s/s0 - (2/pi)arccos(1/eta) ] = K0
    r = (2/np.pi)*np.arccos(1/eta)
    return s0*( K0/(s0*np.sqrt(np.pi*a)*np.sqrt(eta)) + r )

def delta_t(sbar, eta):
    r = sbar/s0 - (2/np.pi)*np.arccos(1/eta)
    return s0*a/Ep*( 4*np.sqrt(eta**2-1)*r + (8/np.pi)*np.log(eta) )

def f(eta):
    return delta_t(sbar_from_ktip(eta), eta) - dc

# find eta where delta_t = dc
etas=np.linspace(1.0001, 5, 2000)
vals=[f(e) for e in etas]
# sign change
br=[]
for i in range(len(etas)-1):
    if vals[i]*vals[i+1]<0: br.append((etas[i],etas[i+1]))
print('brackets eta:',br)
eta=brentq(f, br[0][0], br[0][1], xtol=1e-12)
sbar=sbar_from_ktip(eta)
L=(eta-1)*a
print(f'eta={eta:.8f}  L={L*1e3:.6f} mm  Sigma={sbar/1e6:.6f} MPa')
print(f'  check delta_t={delta_t(sbar,eta)*1e6:.6f} um (dc=9)')
Kt=s0*np.sqrt(np.pi*a)*np.sqrt(eta)*(sbar/s0-(2/np.pi)*np.arccos(1/eta))
print(f'  check K_tip={Kt:.6e} (K0={K0:.6e})')
