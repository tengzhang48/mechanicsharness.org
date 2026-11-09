"""STRONG validation of the field-building + Cartesian FD J machinery:
feed the LINEAR-ELASTIC mode-I field through the SAME code path
(polar -> Cartesian rotation, FD gradients, traction) used for the plastic
field in final.py, and check J = K_I^2 (1-nu^2)/E exactly.
"""
import numpy as np

E = 200e9
nu = 0.3
eps0 = 1.0
sig0 = 1.0          # work in normalized units; exact J is then K^2(1-nu^2)/E
KI = 1.0

# linear-elastic mode-I fields in polar form, r in metres.
def polar_elastic(r, t):
    k = KI / np.sqrt(2 * np.pi * r)
    sr = k * (5.0 / 4.0) * np.cos(t / 2.0) - k * (1.0 / 4.0) * np.cos(3 * t / 2.0)
    st = k * (3.0 / 4.0) * np.cos(t / 2.0) + k * (1.0 / 4.0) * np.cos(3 * t / 2.0)
    srt = k * (1.0 / 4.0) * np.sin(t / 2.0) + k * (1.0 / 4.0) * np.sin(3 * t / 2.0)
    G = E / (2 * (1 + nu))
    kappa = 3 - 4 * nu
    pre = KI / (2 * G) * np.sqrt(r / (2 * np.pi))
    ur = pre * (kappa - np.cos(t)) * np.cos(t / 2.0)
    ut = -pre * (kappa - np.cos(t)) * np.sin(t / 2.0)
    return sr, st, srt, ur, ut


def fold(t):
    tt = np.where(t < 0, t + 2 * np.pi, t)
    return np.where(tt > np.pi, 2 * np.pi - tt, tt)


def sigma_cart(x, y):
    r = np.hypot(x, y); t = np.arctan2(y, x)
    sr, st, srt, ur, ut = polar_elastic(r, t)
    ct = x / r; sy = y / r
    sxx = sr * ct ** 2 + st * sy ** 2 - 2 * srt * ct * sy
    syy = sr * sy ** 2 + st * ct ** 2 + 2 * srt * ct * sy
    sxy = (sr - st) * ct * sy + srt * (ct ** 2 - sy ** 2)
    return sxx, syy, sxy


def u_cart(x, y):
    r = np.hypot(x, y); t = np.arctan2(y, x)
    sr, st, srt, ur, ut = polar_elastic(r, t)
    ct = x / r; sy = y / r
    return ur * ct - ut * sy, ur * sy + ut * ct


def J_numeric(R, nth=80001, h=None):
    if h is None:
        h = 1e-9 * R
    th = np.linspace(-np.pi, np.pi, nth)
    x = R * np.cos(th); y = R * np.sin(th)
    up = u_cart(x + h, y); um = u_cart(x - h, y)
    dux_dx = (up[0] - um[0]) / (2 * h)
    duy_dx = (up[1] - um[1]) / (2 * h)
    sxx, syy, sxy = sigma_cart(x, y)
    nx = np.cos(th); ny = np.sin(th)
    tx = sxx * nx + sxy * ny
    ty = sxy * nx + syy * ny
    # plane-strain energy density
    W = ((1 + nu) / (2 * E)) * ((1 - nu) * (sxx ** 2 + syy ** 2)
                                - 2 * nu * sxx * syy + 2 * sxy ** 2)
    integ = (W * np.cos(th) - (tx * dux_dx + ty * duy_dx)) * R
    return np.trapezoid(integ, th)


for R in [1e-4, 1e-3, 1e-2]:
    Jn = J_numeric(R)
    print(f"R={R:.0e}: J_numeric={Jn:.8e}  K^2(1-nu^2)/E={KI**2*(1-nu**2)/E:.8e}"
          f"  ratio={Jn/(KI**2*(1-nu**2)/E):.8f}")