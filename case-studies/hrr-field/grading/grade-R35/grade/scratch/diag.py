import numpy as np
from verify_elastic2 import sigma_cart, u_cart, KI, E, nu, polar_elastic

# Direct Cartesian stress (verify_elastic.py form)
def Kfield(r, t):
    k = KI / np.sqrt(2 * np.pi * r)
    sxx = k * np.cos(t / 2) * (1 - np.sin(t / 2) * np.sin(3 * t / 2))
    syy = k * np.cos(t / 2) * (1 + np.sin(t / 2) * np.sin(3 * t / 2))
    sxy = k * np.sin(t / 2) * np.cos(t / 2) * np.cos(3 * t / 2)
    return sxx, syy, sxy

for (x, y) in [(1e-3, 0.0), (1e-3, 5e-4), (0.0, 1e-3), (5e-4, 1e-3)]:
    r = np.hypot(x, y); t = np.arctan2(y, x)
    a = sigma_cart(np.array([x]), np.array([y]))
    b = Kfield(r, t)
    print(f"({x:.0e},{y:.0e}): rotated=({a[0][0]:.5e},{a[1][0]:.5e},{a[2][0]:.5e})  direct=({b[0]:.5e},{b[1]:.5e},{b[2]:.5e})")

# Displacement gradient check: compare FD of u_cart to analytic eps at theta=0.
# At theta=0 (y=0, x>0): ux = ur(x,0)*(1) - ut*0 = ur; ur = pre*(kappa-cos0)*cos0
# = pre*(kappa-1); pre = KI/(2G) sqrt(r/2pi). So ux = (KI/(2G))sqrt(x/2pi)(kappa-1)
from numpy import sqrt, pi
x = 1e-3; h = 1e-9
def ux_analytic(xx):
    pre = KI / (2 * (E / (2 * (1 + nu)))) * np.sqrt(xx / (2 * pi))
    return pre * (3 - 4 * nu - 1)
up = u_cart(np.array([x + h]), np.array([0.0]))[0][0]
um = u_cart(np.array([x - h]), np.array([0.0]))[0][0]
print("du_x/dx (FD from u_cart):", (up - um) / (2 * h))
print("du_x/dx (analytic):", (ux_analytic(x + h) - ux_analytic(x - h)) / (2 * h))
# also ux at theta=0 should equal ur
print("u_x(x,0) from u_cart:", u_cart(np.array([x]), np.array([0.0]))[0][0], " ur analytic:", ux_analytic(x))