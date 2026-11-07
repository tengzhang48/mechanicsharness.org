import numpy as np, march
from march import _P1, _P2, _P3

# Problem statement, dimensionless, S(x,y) with traction s = 1 - (1/pi) PV Int S(x,y) w(y) dy:
#   S(x,y) = 1/(x-y) - 1/(x+y) - 2y(x-y)/(x+y)^3
def S_stmt(x, y):
    return 1/(x-y) - 1/(x+y) - 2*y*(x-y)/(x+y)**3

# What march.py actually implements (d/dy of its primitives with c0=1,c1=0):
def K_impl(x, y):
    return 1/(y-x) + 1/(y+x) - 4*x/(y+x)**2 + 4*x**2/(y+x)**3

# Simplify S_stmt:  -2y(x-y)/(x+y)^3 = -2y(x+y-2y)/(x+y)^3 = -2y/(x+y)^2 + 4y^2/(x+y)^3
def S_stmt2(x, y):
    return 1/(x-y) - 1/(x+y) - 2*y/(x+y)**2 + 4*y**2/(x+y)**3

# K_impl: 1/(y-x) = -1/(x-y); so K_impl = -1/(x-y) + 1/(x+y) - 4x/(x+y)^2 + 4x^2/(x+y)^3
# These are clearly DIFFERENT kernels (signs & which variable multiplies).
for (x, y) in [(1.2, 0.3), (0.3, 1.2), (0.7, 0.2)]:
    print("x=%.2f y=%.2f  S_stmt=%.6f  S_stmt2=%.6f  K_impl=%.6f"
          % (x, y, S_stmt(x, y), S_stmt2(x, y), K_impl(x, y)))