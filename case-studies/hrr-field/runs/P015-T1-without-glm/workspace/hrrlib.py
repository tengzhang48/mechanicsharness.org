"""
HRR plane-strain crack-tip field library, n = 6.5.
Airy: Phi = r^(2-p) phi(theta):
  sig_rr = (2-p)phi + phi'', sig_tt = (2-p)(1-p)phi, sig_rt = -(1-p)phi'
plane strain: sig_zz=(sig_rr+sig_tt)/2, s_rr=(sig_rr-sig_tt)/2=-s_tt, s_rt=sig_rt,
  sigma_e = sqrt(3(s_rr^2+s_rt^2)).
eps ~ r^(-pn) with p(n+1)=1; u_r = r^p a, u_th = r^p b give
  e~_rr = p a, e~_tt = a+b', 2 e~_rt = a' + (p-1)b, and e~_tt = -e~_rr, so
  ODE: (e~_rr)'' + (1-p^2) e~_rr - 2p (e~_rt)' = 0,
  e~_rr = (3/2) se^(n-1) s_rr, e~_rt = (3/2) se^(n-1) s_rt, se = sig~_e(theta).
BCs: phi'(0)=0, phi'''(0)=0 (mode-I symmetry); phi(pi)=phi'(pi)=0 (traction free).
J = alpha e0 sigma0 I_n (B/sigma0)^(n+1); sigma = sigma0 (J/(alpha e0 sigma0 I_n r))^p sig~
"""
import numpy as np
from scipy.integrate import solve_bvp

PI = np.pi
n_t = 6.5
alpha, E, s0 = 0.02, 200e9, 400e6
e0 = s0 / E

def _derive_F():
    import sympy as sp
    t, ns, ps = sp.symbols('t ns ps', positive=True)
    phi = sp.Function('PHI')(t)
    y1, y2, y3, y4, y5 = sp.symbols('y1 y2 y3 y4 y5')
    s1 = (sp.diff(phi, t, 2) + ps * (2 - ps) * phi) / 2
    s2 = -(1 - ps) * sp.diff(phi, t, 1)
    se = sp.sqrt(3 * (s1 ** 2 + s2 ** 2))
    q1 = se ** (ns - 1) * s1
    q2 = se ** (ns - 1) * s2
    expr = sp.diff(q1, t, 2) + (1 - ps ** 2) * q1 - 2 * ps * sp.diff(q2, t, 1)
    sub = {sp.Derivative(phi, (t, 1)): y2, sp.Derivative(phi, (t, 2)): y3,
           sp.Derivative(phi, (t, 3)): y4, sp.Derivative(phi, (t, 4)): y5, phi: y1}
    w = sp.solve(sp.Eq(expr.subs(sub), 0), y5)[0]
    return sp.lambdify((y1, y2, y3, y4, ns, ps), w, 'numpy')

_FSYM = _derive_F()

def F_phi4(y1, y2, y3, y4, n, p):
    """Fourth derivative phi'''' from the compatibility ODE (sympy-derived)."""
    y1, y2, y3, y4 = float(y1), float(y2), float(y3), float(y4)
    if n == 1.0:
        s1 = (y3 + p * (2.0 - p) * y1) / 2.0
        s2p = -(1.0 - p) * y3
        return -2.0 * (p * (2.0 - p) * y3 / 2.0 + (1.0 - p * p) * s1 - 2.0 * p * s2p)
    return float(_FSYM(y1, y2, y3, y4, n, p))

def make_bvp(nv):
    def ode(th, y, pp):
        p = pp[0]
        return np.vstack([y[1], y[2], y[3],
                          [_FSYM(y[0], y[1], y[2], y[3], nv, p)]]) * 1.0
    def bc(ya, yb, pp):
        return np.array([ya[1], ya[3], yb[0], yb[1], ya[0] - 1.0])
    return ode, bc
