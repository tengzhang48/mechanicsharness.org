"""Algebraic traction route: calibrate the large-p tail relation on the b1 column.
T(p) = M(p,1) * Phi(p); for b1 (X^-1/2 e^-X): Phi = sqrt(pi) (1+ip)^{-1/2}.
Known: A32 = -cc/2 (xx), A12 = cc/2 (xx), and for xz-column: A12 = -i cc nu.
Tail model: M Phi = C32 * p^{1/2} + C12 * p^{-1/2} + O(p^{-3/2}); find C32, C12
numerically at large real p, then determine A = C * const (calibration)."""
import numpy as np
from scipy.special import gamma as G
from m_stiffness import M_closed

NU = 0.3; CC = -1.0/(2*(1-NU))

ps = np.array([2e3, 4e3, 8e3, 1.6e4, 3.2e4])
def tails(entry_fn, ps):
    """returns fitted C32, C12 of T(p) = C32 p^1/2 + C12 p^-1/2 (+c0 p^-3/2)"""
    T = np.array([entry_fn(p) for p in ps])
    Mfit = np.stack([ps**0.5, ps**-0.5, ps**-1.5], axis=1)
    c = np.linalg.lstsq(Mfit, T, rcond=None)[0]
    return c[0], c[1]

# b1 column x (v = e_x): Phi = sqrt(pi)(1+ip)^{-1/2}
def Txx_b1x(p):
    return M_closed(p, 1.0)[0, 0]*np.sqrt(np.pi)*(1+1j*p)**-0.5
def Tzz_b1x(p):
    return M_closed(p, 1.0)[2, 2]*np.sqrt(np.pi)*(1+1j*p)**-0.5
def Tzz_b1z(p):
    return M_closed(p, 1.0)[2, 2]*np.sqrt(np.pi)*(1+1j*p)**-0.5   # same Phi (unit edge)
def Txz_b1z(p):
    return M_closed(p, 1.0)[0, 2]*np.sqrt(np.pi)*(1+1j*p)**-0.5

C32xx, C12xx = tails(Txx_b1x, ps)
C32zz, C12zz = tails(Tzz_b1x, ps)
_, C12xz = tails(Txz_b1z, ps)
A32xx_ex, A12xx_ex = -CC/2, CC/2
A32zz_ex, A12zz_ex = -CC*(1-NU)/2, CC*(1-NU)/2
A12xz_ex = -1j*CC*NU
print("calibration ratios (should be constant across entries):")
print(f"  xx: A32/C32 = {A32xx_ex/C32xx:.8f}   A12/C12 = {A12xx_ex/C12xx:.8f}")
print(f"  zz: A32/C32 = {A32zz_ex/C32zz:.8f}   A12/C12 = {A12zz_ex/C12zz:.8f}")
print(f"  xz: A12/C12 = {A12xz_ex/C12xz:.8f}")
print("absolute:")
print(f"  C32xx={C32xx:.6f} A32xx={A32xx_ex:.6f} | C12xx={C12xx:.6f} A12xx={A12xx_ex:.6f}")
print(f"  C12xz={C12xz:.6f} A12xz={A12xz_ex:.6f}")
