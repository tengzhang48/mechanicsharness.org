import numpy as np
from kernels_x import BasisFun, apply_ahead_matrix, apply_kernels_matrix

# ahead traction of the raw (unscaled) m=0 minus basis function: expect
# T(x) ~ A32/|x|^{3/2} + ... with A32 = 1/(sqrt(2 pi) D1) = 0.35714 (edge coeff 1)
bf = BasisFun(-0.5, 0)
xs = np.array([3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1])
Tv = np.array([apply_ahead_matrix([bf], -x, ['yy'])['yy'][0] for x in xs])
print("ahead T of X^-1/2 e^-X:", Tv.real)
print("implied A32 = T*x^1.5 :", (Tv*xs**1.5).real)
# on-crack traction of same (should be far from 0)
print("on-crack T at 0.1, 1.0:", apply_kernels_matrix([bf], 0.1, ['yy'])['yy'][0],
      apply_kernels_matrix([bf], 1.0, ['yy'])['yy'][0])
