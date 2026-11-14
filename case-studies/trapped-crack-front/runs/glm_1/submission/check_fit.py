"""Check ahead-evaluation accuracy at small |x| and the sqrt(x)-extrapolation fit:
T_reg(x) x^{1/2} = A12 + A0 u + A1 u^2 + ... (u = sqrt(x)).
Test with the m=0 plus function whose A12 can be cross-checked via cut quadrature.
"""
import numpy as np
from kernels_x import BasisFun, apply_ahead_matrix

bf = BasisFun(0.5, 0)
xs = np.array([1e-6, 4e-6, 1.6e-5, 6.4e-5, 2.56e-4, 1.024e-3, 4.1e-3, 1.6e-2])
for entry in ['xx']:
    Tv = np.array([apply_ahead_matrix([bf], -x, [entry])[entry][0] for x in xs])
    u = np.sqrt(xs)
    fx = Tv*u           # = A12 + A0 u + ...
    print("f(u) = T*sqrt(x):", np.array_str(fx.real, precision=8))
    for deg in [1, 2, 3, 4]:
        c = np.polyfit(u, fx.real, deg)
        print(f"  deg {deg}: A12 = {np.polyval(c, 0.0):.8f}")
    # also T*x^{3/2} -> A32 (expect ~ 0)
    g = Tv*u**3
    for deg in [2, 3, 4]:
        c = np.polyfit(u, g.real, deg)
        print(f"  A32 deg {deg}: {np.polyval(c, 0.0):.6e}")
