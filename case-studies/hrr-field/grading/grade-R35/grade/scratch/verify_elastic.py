"""Validate the Cartesian finite-difference J machinery on the elastic field,
where the exact answer is J = K_I^2/E (plane stress)."""
import numpy as np

alpha = 0.02
sig0 = 400e6
E = 200e9
eps0 = sig0 / E
nu = 0.3

# Elastic mode-I angular fields (plane stress), s = 3/2.
# phi = K_I/sqrt(2 pi) r^{3/2}[cos(t/2)+(1/3)cos(3t/2)]; stresses sigma = ...
# Standard mode-I: K_I amplitude. Use the KNOWN closed forms.
def Kfield(KI, r, t):
    k = KI / np.sqrt(2 * np.pi * r)
    ct = np.cos(t / 2); s2 = np.sin(t) / 2.0
    sxx = k * np.cos(t / 2) * (1 - np.sin(t / 2) * np.sin(3 * t / 2))
    syy = k * np.cos(t / 2) * (1 + np.sin(t / 2) * np.sin(3 * t / 2))
    sxy = k * np.sin(t / 2) * np.cos(t / 2) * np.cos(3 * t / 2)
    sxy = k * np.sin(t / 2) * np.cos(t / 2) * np.cos(3 * t / 2)
    return sxx, syy, sxy

# displacements (plane strain) with G=E/(2(1+nu))
def ufield(KI, r, t, E, nu):
    G = E / (2 * (1 + nu))
    kappa = 3 - 4 * nu
    pre = KI / (2 * G) * np.sqrt(r / (2 * np.pi))
    ux = pre * (kappa - np.cos(t)) * np.cos(t / 2)
    uy = pre * (kappa - np.cos(t)) * np.sin(t / 2)
    return ux, uy

def J_cart(KI, R, E, nu, nth=20001, h=1e-8):
    th = np.linspace(-np.pi, np.pi, nth)
    x = R * np.cos(th); y = R * np.sin(th)
    t = np.arctan2(y, x)
    # displacement gradients in x
    def u(xx, yy):
        tt = np.arctan2(yy, xx)
        return ufield(KI, np.hypot(xx, yy), tt, E, nu)
    up = u(x + h, y); um = u(x - h, y)
    dux_dx = (up[0] - um[0]) / (2 * h)
    duy_dx = (up[1] - um[1]) / (2 * h)
    sxx, syy, sxy = Kfield(KI, R + 0 * th, th)
    nx = np.cos(th); ny = np.sin(th)
    tx = sxx * nx + sxy * ny
    ty = sxy * nx + syy * ny
    # plane-strain W = (1+nu)/(2E)[(1-nu)(sxx^2+syy^2) - 2 nu sxx syy + 2 sxy^2]
    W = ((1 + nu) / (2 * E)) * ((1 - nu) * (sxx ** 2 + syy ** 2)
                                - 2 * nu * sxx * syy + 2 * sxy ** 2)
    integ = (W * np.cos(th) - (tx * dux_dx + ty * duy_dx)) * R
    return np.trapezoid(integ, th)

for R in [1e-3, 1e-2, 1e-1]:
    Jn = J_cart(1e5, R, E, nu)
    print(f"R={R:.0e}: J_numeric={Jn:.6e}  K^2/E={1e5**2/E:.6e}  K^2(1-nu^2)/E={1e5**2*(1-nu**2)/E:.6e}  ratio_ps={Jn/(1e5**2/E):.6f}")