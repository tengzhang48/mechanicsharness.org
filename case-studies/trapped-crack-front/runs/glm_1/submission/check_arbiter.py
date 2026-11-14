"""Arbiter for T(b = X^{3/2} e^{-X}) ahead (x<0), y-sector:
(1) X-domain kernel integral (kernels_x.apply_ahead_matrix)
(2) p-domain QAWO real-line
(3) lower-cut contour (validated constant -i, no finite parts needed here)
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma as G
from kernels_x import BasisFun, apply_ahead_matrix
from m_stiffness import M_closed

NU = 0.3; MU = 1.0; cc = -MU/(2*(1-NU)); gam = G(2.5)

def T_contour(x):
    # F = -(i/2pi) int (fR - fL) e^{-(1+s)|x|} ds ; f_pm = cc*(+- i rho) * Phi_low
    # fR - fL = 2 i cc rho * Phi_low,  Phi_low = gam*(2+s)^{-5/2}
    f = lambda s: 2j*cc*np.sqrt(s*(s+2.0))*gam*(2.0+s)**-2.5
    I = quad(lambda s: f(s).real, 0, 200, limit=400)[0] + \
        1j*quad(lambda s: f(s).imag, 0, 200, limit=400)[0]
    return -1j/(2*np.pi)*I*np.exp(0)  # note e^{-(1+s)|x|} inside integrand below

def T_contour2(x):
    f = lambda s: (2j*cc*np.sqrt(s*(s+2.0))*gam*(2.0+s)**-2.5)*np.exp(-(1.0+s)*abs(x))
    I = quad(lambda s: f(s).real, 0, 200, limit=400)[0] + \
        1j*quad(lambda s: f(s).imag, 0, 200, limit=400)[0]
    return -1j/(2*np.pi)*I

def T_qawo(x):
    f = lambda p: (M_closed(p, 1.0)[1, 1]*gam*(1+1j*p)**-2.5)
    Ic = quad(lambda p: f(p).real, 0, np.inf, weight='cos', wvar=x, limlst=400)[0]
    Is = quad(lambda p: f(p).imag, 0, np.inf, weight='sin', wvar=x, limlst=400)[0]
    return 2*(Ic - Is)/(2*np.pi)

bf = BasisFun(1.5, 0)
for x in [0.1, 0.5, 2.0]:
    v1 = apply_ahead_matrix([bf], -x, ['yy'])['yy'][0]
    v2 = T_qawo(x)
    v3 = T_contour2(x)
    print(f"x=-{x:5.3f}: X-domain={v1.real:+.8f}  QAWO={v2.real:+.8f}  contour={v3.real:+.8f}")
