import numpy as np
from numpy.polynomial.legendre import leggauss

def J0_invsqrt(x):
    d = np.sqrt(1.0 - x)
    return np.log((1.0 - d) / (1.0 + d)) / d

def quad(lo, hi, func, n=3000):
    z, wz = leggauss(n)
    m = (lo + hi) / 2; r = (hi - lo) / 2
    return r * np.sum(wz * func(m + r * z))

for x in [0.05, 0.2, 5.0e-1, 0.8, 0.95]:
    d = np.sqrt(1.0 - x)
    # 2/(t^2-d^2) = (1/d)[1/(t-d) - 1/(t+d)]
    # regular part: -(1/d)/(t+d)
    pv = -(1.0 / d) * (quad(0, 1, lambda u: 1.0 / (u + d)))
    # PV int_0^1 dt/(t-d) = ln(1-d) - ln(d)
    pv += (1.0 / d) * (np.log(1 - d) - np.log(d))
    exact = J0_invsqrt(x)
    print(f"x={x}: closed={exact:.14f} brute={pv:.14f} diff={exact-pv:.3e}")
