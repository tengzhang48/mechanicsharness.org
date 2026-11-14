"""(a) Validate the plus-family m=0 function (nu=3/2) ahead values: Hadamard
machinery vs lower-cut quadrature (Abel).  (b) Residual profile of the shear LS.
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma as G
from kernels_x import BasisFun, apply_ahead_matrix, apply_kernels_matrix

CC = -1/(2*0.7); NU = 0.3

def Mentry(entry, p, r):
    if entry == 'xx': return CC*((1-NU)*r + NU*p*p/r)
    if entry == 'zz': return CC*((1-NU)*r + NU/r)
    if entry == 'xz': return CC*(NU*p/r)

def Phat_Qhat(entry, s):
    if entry == 'xx': return CC*NU*np.ones_like(s), CC*np.ones_like(s)
    if entry == 'zz': return -CC*NU*np.ones_like(s), CC*(1-NU)*np.ones_like(s)
    if entry == 'xz': return 1j*CC*NU*(1.0+s), np.zeros_like(s)

def T_cut_quad_gen(entry, x, alpha, m=0):
    """ahead Abel traction of e^{-X} X^alpha L_m(X) via lower cut:
    Phi_low(p) at p=-i(1+s): sum_k c_k Gamma(alpha+k+1) (2+s)^{-(alpha+k+1)}."""
    from math import comb, factorial
    c = np.array([(-1)**k*comb(m, k)/factorial(k) for k in range(m+1)])
    expo = alpha + np.arange(m+1) + 1
    def Phi(s):
        return sum(c[k]*G(alpha+k+1)*(2.0+s)**(-(alpha+k+1)) for k in range(m+1))
    def g(s):
        P, Q = Phat_Qhat(entry, s)
        rho = np.sqrt(s*(s+2.0))
        fRmL = -2*((1j/rho)*P + 1j*rho*Q)*Phi(s)
        return fRmL*np.exp(-(1.0+s)*x)
    val = quad(lambda s: np.real(g(s)), 0, 400, limit=500)[0] \
        + 1j*quad(lambda s: np.imag(g(s)), 0, 400, limit=500)[0]
    return -1j/(2*np.pi)*val

bf = BasisFun(0.5, 0)
print("ahead T of X^{1/2} e^{-X} (nu=3/2): Hadamard vs cut-quad:")
for entry in ['xx', 'zz', 'xz']:
    for x in [0.1, 0.5, 2.0]:
        v1 = apply_ahead_matrix([bf], -x, [entry])[entry][0]
        v2 = T_cut_quad_gen(entry, x, 0.5)
        print(f"  {entry} x={x:4.1f}: hadamard={v1:+.8f}  cut={v2:+.8f}  diff={v1-v2:+.2e}")
print("on-crack check at X=0.5, 2.0 (hadamard only, for reference):")
for entry in ['xx', 'zz']:
    print(" ", entry, apply_kernels_matrix([bf], 0.5, [entry])[entry][0],
          apply_kernels_matrix([bf], 2.0, [entry])[entry][0])
