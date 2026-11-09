"""
Two unequal strips joined along a weak plane; opening couples M at the free ends
(no end forces); starter cut 20 mm; newly separated faces carry softening
connectors t(delta) = sigma_b (1 - delta/delta_c), 0<delta<delta_c, inactive
beyond delta_c.  Front advance consumes Gamma0 (J/m^2 of joining-plane area).

Model (Euler-Bernoulli; Bao & Suo 1992 splitting-beam setup):
 - s runs from the advancing front (s=0) toward the starter-cut root (s=L).
 - Bonded part ahead: long, straight -> arms clamped at the front:
   w_i(0)=0, w_i'(0)=0  =>  delta(0)=0, delta'(0)=0.
 - Behind the root (starter cut, no connectors, no end forces): each arm in
   pure bending, moment M, zero shear -> delta''(L)=C_M*M, delta'''(L)=0,
   C_M = 1/D1 + 1/D2, D_i = E_i I_i.
 - Arms: D_i w_i'''' = -b t(delta)  =>  delta'''' = mu^4 (delta - delta_c),
   mu^4 = b C_M sigma_b / delta_c.   (softening -> negative stiffness)
 - Growth condition (exact Rice J-integral; Bao&Suo eqs (15),(22),(51)):
       J_far = M^2 C_M/(2b) = Gamma0 + int_0^{delta_t} t dd,
   delta_t = delta(L).  Event: first connector reaches delta_c  <=>  delta_t=delta_c
   (delta(s) monotone, max at s=L).
"""
import numpy as np
from scipy.optimize import brentq

# ---------------- constants (SI) ----------------
b   = 25.0e-3
E1, E2 = 110e9, 190e9
h1, h2 = 0.60e-3, 0.90e-3
sig  = 0.20e6
dc   = 0.20e-3
G0   = 1.40
ell0 = 20.0e-3

D1 = E1 * b * h1**3 / 12
D2 = E2 * b * h2**3 / 12
CM = 1/D1 + 1/D2
mu4 = b * CM * sig / dc
mu  = mu4**0.25

def basis(s):
    """phi_k, phi', phi'', phi''' (k: cosh, sinh, cos, sin) of y''''=mu^4 y."""
    x = mu*np.asarray(s, float)
    ch, sh, co, si = np.cosh(x), np.sinh(x), np.cos(x), np.sin(x)
    m = mu
    P   = np.array([ch, sh, co, si])
    dP  = m*np.array([sh, ch, -si, co])
    d2P = m**2*np.array([ch, sh, -co, -si])
    d3P = m**3*np.array([sh, ch, si, -co])
    return P, dP, d2P, d3P

def solution(L, M, s):
    """delta(s); BCs delta(0)=0, delta'(0)=0, delta''(L)=CM*M, delta'''(L)=0."""
    P, dP, d2P, d3P = [np.asarray(v).reshape(4) for v in basis(L)]
    # delta = dc + c0 cosh + c1 sinh + c2 cos + c3 sin
    # delta(0)=0:  dc + c0 + c2 = 0
    # delta'(0)=0: c1 + c3 = 0
    # rows at L for unknowns (c0, c1) with c2 = -dc - c0, c3 = -c1:
    #   c0(A2[0]-A2[2]) + c1(A2[1]-A2[3]) = CM*M + dc*A2[2]
    #   c0(A3[0]-A3[2]) + c1(A3[1]-A3[3]) = dc*A3[2]
    a11, a12 = d2P[0]-d2P[2], d2P[1]-d2P[3]
    a21, a22 = d3P[0]-d3P[2], d3P[1]-d3P[3]
    r1 = CM*M + dc*d2P[2]
    r2 = dc*d3P[2]
    det = a11*a22 - a12*a21
    c0 = (r1*a22 - a12*r2)/det
    c1 = (a11*r2 - a21*r1)/det
    c = np.array([c0, c1, -dc-c0, -c1])
    Pb, *_ = basis(s)
    d = dc + c @ Pb
    return (d if np.ndim(d) else np.atleast_1d(d)), c

def delta_end(L, M):
    return float(np.asarray(solution(L, M, [L])[0]).reshape(-1)[0])



def bridging_work(dt):
    return sig*dt*(1 - 0.5*dt/dc)
Minit = np.sqrt(2*b*G0/CM)
Mstar = np.sqrt(2*b*(G0 + sig*dc/2)/CM)
