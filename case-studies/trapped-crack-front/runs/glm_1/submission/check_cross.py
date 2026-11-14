"""Check: is b = X^{-1/2} e^{-X} exactly traction-free on the crack (y-sector)?
Contour cut analysis says fR-fL = 0 identically (cc*r * sqrt(pi)(1+ip)^-1/2:
M-jump and Phi-jump cancel).  X-domain Hadamard machinery gives independent verdict.
Also check the ahead traction (should be nonzero)."""
import numpy as np
from kernels_x import BasisFun, apply_kernels_matrix, apply_ahead_matrix

bf = BasisFun(-0.5, 0)
print("on-crack T_yy(X) for b = X^-1/2 e^-X:")
for X in [0.02, 0.05, 0.2, 1.0, 3.0, 10.0]:
    v = apply_kernels_matrix([bf], X, ['yy'])['yy'][0]
    print(f"  X={X:8.3f}:  T = {v.real:+.3e} {v.imag:+.3e}j")
print("ahead T_yy(x) for x<0:")
for x in [0.02, 0.1, 0.5, 2.0]:
    v = apply_ahead_matrix([bf], -x, ['yy'])['yy'][0]
    print(f"  x=-{x:6.3f}:  T = {v.real:+.6e} {v.imag:+.3e}j")
# also shear sector pure functions
bfz = BasisFun(-0.5, 0)
print("\nshear: on-crack T of b_x = X^-1/2 e^-X (xx and xz entries):")
for X in [0.1, 1.0]:
    d = apply_kernels_matrix([bf], X, ['xx', 'xz'])
    print(f"  X={X}: Txx = {d['xx'][0]:+.6e}  Txz = {d['xz'][0]:+.6e}")
