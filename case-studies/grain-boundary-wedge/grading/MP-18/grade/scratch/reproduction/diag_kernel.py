import numpy as np, march
from march import _P1, _P2, _P3

zc = 1.2345
h = 1e-6
for yv in [2.13, 0.41, 0.91]:
    def F(y): return _P1(y, zc) + _P2(y, zc) + _P3(y, zc)
    dF = (F(yv + h) - F(yv - h)) / (2 * h)
    print("y=%.2f  d/dy(P1+P2+P3)=%+.8f" % (yv, dF))

# Also derive analytically:
# d/dy log|y-z| = 1/(y-z)
# d/dy log|y+z| = 1/(y+z)
# d/dy [4z/(y+z) - 2z^2/(y+z)^2] = -4z/(y+z)^2 + 4z^2/(y+z)^3
# so latent kernel implied by primitives:
def K_implied(z, y):
    return 1/(y-z) + 1/(y+z) - 4*z/(y+z)**2 + 4*z**2/(y+z)**3

def K_stmt(z, y):   # problem statement
    return 1/(z-y) - 1/(z+y) - 2*y*(z-y)/(z+y)**3

def K_march_doc(z, y):
    return 1/(y-z) + 1/(y+z) - 2*z*(y-z)/(y+z)**3

for (z, y) in [(1.2, 0.3), (1.2, 2.0), (0.5, 0.5 + 0.1)]:
    print("z=%.2f y=%.2f  implied=%.6f  stmt=%.6f  march_doc=%.6f"
          % (z, y, K_implied(z, y), K_stmt(z, y), K_march_doc(z, y)))