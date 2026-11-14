"""Test the algebraic tail route on phi0+ (known A12[phi0+,xx] = 0.357138 from the
validated sqrt(x)-fit, which in turn matches the exact -A12[b1] = -cc/2 identity),
and A12[phi0+,xz] = 0 (validated fit, residual 1e-7)."""
import numpy as np
from scipy.special import gamma as G
from m_stiffness import M_closed

CAL12 = (1-1j)/np.sqrt(2*np.pi)
ps = np.array([2e3, 4e3, 8e3, 1.6e4, 3.2e4])
Mfit = np.stack([ps**0.5, ps**-0.5, ps**-1.5, ps**-2.5], axis=1)

def A12_tail(coef_fun):
    T = np.array([coef_fun(p) for p in ps])
    c = np.linalg.lstsq(Mfit, T, rcond=None)[0]
    return CAL12*c[1]

gam32 = G(1.5)
# phi0+ in x-component: Phi = Gamma(3/2)(1+ip)^{-3/2}
print("A12[phi0+,xx] algebraic-tail:", A12_tail(lambda p: M_closed(p,1.0)[0,0]*gam32*(1+1j*p)**-1.5), " expect +0.357138")
print("A12[phi0+,zz] algebraic-tail:", A12_tail(lambda p: M_closed(p,1.0)[2,2]*gam32*(1+1j*p)**-1.5), " expect +0.249997")
print("A12[phi0+,xz] algebraic-tail:", A12_tail(lambda p: M_closed(p,1.0)[0,2]*gam32*(1+1j*p)**-1.5), " expect 0")
