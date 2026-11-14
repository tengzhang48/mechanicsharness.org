"""Compute A12[phi0+] for all kernel entries (phi0+ = X^{1/2} e^{-X}, unit sqrt-edge)
and compare with -A12[b1] per entry:
  b1 (X^{-1/2}e^{-X}):  A12^b1_xx = cc/2, A12^b1_zz = cc(1-nu)/2, A12^b1_xz = -i*cc*nu
"""
import numpy as np
from kernels_x import BasisFun, apply_ahead_matrix

CC = -1/(2*0.7); NU = 0.3
bf = BasisFun(0.5, 0)
xs = np.array([1e-6, 4e-6, 1.6e-5, 6.4e-5, 2.56e-4, 1.024e-3, 4.1e-3, 1.6e-2])
u = np.sqrt(xs)
M12 = np.stack([np.ones_like(u), u**2, u**4, u**6], axis=1)
print("A12 of phi0+ (unit sqrt-edge X^{1/2} e^{-X}):")
for entry in ['xx', 'zz', 'xz']:
    Tv = np.array([apply_ahead_matrix([bf], -x, [entry])[entry][0] for x in xs])
    s = np.linalg.lstsq(M12, Tv*u**2, rcond=None)[0]
    A12 = s[0]
    pred = {'xx': -CC/2, 'zz': -CC*(1-NU)/2, 'xz': +1j*CC*NU}[entry]
    print(f"  {entry}: A12 = {A12:+.8f}   -A12^b1 = {pred:+.8f}   ratio={A12/pred:.8f}")
