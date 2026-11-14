"""Decisive check of the X-domain Hadamard kernel application:
for b = X^{3/2} e^{-X}, M_yy*Phat ~ p^{-3/2} decays, so
T(X) = (1/2pi) int M_yy(p,1) Gamma(5/2)(1+ip)^{-5/2} e^{ipX} dp
converges and is computed by direct oscillatory quadrature (QAWO)."""
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma as G
from kernels_x import BasisFun, apply_kernels_matrix, apply_ahead_matrix, kernels
from m_stiffness import M_closed

NU = 0.3; MU = 1.0
gam = G(2.5)

def T_pdomain(X):
    f = lambda p: (M_closed(p, 1.0)[1, 1]*gam*(1+1j*p)**-2.5)
    re = lambda p: (f(p)).real
    im = lambda p: (f(p)).imag
    Ic = quad(lambda p: re(p)*np.cos(0), 0, 0)  # placeholder
    # even/odd split: f(-p) = conj(f(p)) (M_yy even, (1+ip)^{-5/2} conj-symmetric)
    # int_{-inf}^{inf} f e^{ipX} dp = 2 int_0^inf [Re f cos(pX) - Im f sin(pX)] dp
    Ic = quad(lambda p: re(p), 0, np.inf, weight='cos', wvar=X, limlst=400)[0]
    Is = quad(lambda p: im(p), 0, np.inf, weight='sin', wvar=X, limlst=400)[0]
    return 2*(Ic - Is)/(2*np.pi)

bf = BasisFun(1.5, 0)
print("on-crack T_yy(X): X-domain Hadamard vs p-domain QAWO")
for X in [0.05, 0.2, 1.0, 3.0]:
    v1 = apply_kernels_matrix([bf], X, ['yy'])['yy'][0]
    v2 = T_pdomain(X)
    print(f"  X={X:6.3f}:  hadamard={v1.real:+.8f}   p-domain={v2.real:+.8f}   diff={v1.real-v2.real:+.2e}")
print("ahead T_yy(x), x<0:")
def T_pdomain_neg(x):
    f = lambda p: (M_closed(p, 1.0)[1, 1]*gam*(1+1j*p)**-2.5)
    Ic = quad(lambda p: f(p).real, 0, np.inf, weight='cos', wvar=x, limlst=400)[0]
    Is = quad(lambda p: f(p).imag, 0, np.inf, weight='sin', wvar=x, limlst=400)[0]
    # int f e^{ipx} dp with x<0: = 2 int_0 [Re f cos - i Im f * (-sin)?] : e^{ipx}, x<0
    # f(-p)=conj f(p): int_{-inf}^{inf} = 2 int_0^inf [Re f cos(px) - Im f sin(px)] dp (same formula, x<0)
    return 2*(Ic - Is)/(2*np.pi)
for x in [0.1, 0.5, 2.0]:
    v1 = apply_ahead_matrix([bf], -x, ['yy'])['yy'][0]
    v2 = T_pdomain_neg(x)
    print(f"  x=-{x:5.3f}:  hadamard={v1.real:+.8f}   p-domain={v2.real:+.8f}   diff={v1.real-v2.real:+.2e}")
