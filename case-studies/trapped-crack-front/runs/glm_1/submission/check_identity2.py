"""A12 of phi0+ (X^{1/2} e^{-X}) for the xz kernel with the CORRECT fit model,
vs -A12^{b1,xz} = -(-i cc nu) = -0.214285...i ?   A12^{b1,xz} = -i*cc*nu = +0.2143i.
Identity A12[phi0+] = -A12[b1] would require A12[phi0+,xz] = -0.2143i.
Also re-check xx: expect +0.3571429.
"""
import numpy as np
from kernels_x import BasisFun, apply_ahead_matrix
from first_order import ahead_coeffs_reg, b1_A32_A12

bf = BasisFun(0.5, 0)
for entry in ['xx', 'zz', 'xz']:
    A32, A12, r = ahead_coeffs_reg([bf], np.ones(1), entry)
    _, A12b1 = b1_A32_A12(entry)
    print(f"{entry}: A12[phi0+] = {A12:+.8f}   -A12[b1] = {-A12b1:+.8f}  "
          f"ratio={A12/-A12b1:.6f}  (A32fit={A32:.1e}, fitres={r:.1e})")
