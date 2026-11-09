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

# ---------------- event ----------------
def bridging_work(dt):        # int_0^dt t dd
    return sig*dt*(1 - 0.5*dt/dc)

Mstar = np.sqrt(2*b*(G0 + bridging_work(dc))/CM)
Jstar = Mstar**2 * CM / (2*b)

f = lambda L: delta_end(L, Mstar) - dc
Ls = np.linspace(1e-4, 60e-3, 6000)
vals = np.array([f(L) for L in Ls])
print("sign changes at L(mm):", Ls[np.nonzero(np.diff(np.sign(vals)))[0]]*1e3)
i = np.nonzero(np.diff(np.sign(vals)))[0][0]
Lstar = brentq(f, Ls[i], Ls[i+1], xtol=1e-16, rtol=8.9e-16)
d_mid = float(np.asarray(solution(Lstar, Mstar, [Lstar/2])[0]).reshape(-1)[0])

print(f"D1 = {D1:.6e}, D2 = {D2:.6e} N m^2; C_M = {CM:.8f} 1/(N m^2)")
print(f"mu = {mu:.6f} 1/m (1/mu = {1/mu*1e3:.4f} mm)")
print(f"plateau bridging work sigma_b*delta_c/2 = {bridging_work(dc):.3f} J/m^2")
print(f"J* = {Jstar:.6f} J/m^2")
print(f"M* = {Mstar:.8f} N m")
print(f"L* = {Lstar*1e3:.6f} mm")
print(f"delta(L*/2) = {d_mid*1e6:.6f} um")

# ---------------- checks ----------------
print("\n--- checks ---")
ss = np.linspace(0, Lstar, 5001)
prof, _ = solution(Lstar, Mstar, ss)
print("profile monotone:", bool(np.all(np.diff(prof) > 0)),
      "; max interior (um):", prof.max()*1e6, "; <= delta_c:", prof.max() <= dc+1e-12)
print("bending-only comparison: delta(L/2)_beam = ",
      0.5*CM*Mstar*(Lstar/2)**2*1e6, "um (springs must reduce this)")

# 1) independent ODE integration
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
def rhs_ode(s, y):
    return [y[1], y[2], y[3], mu4*(y[0]-dc)]
def shoot(par):
    so = solve_ivp(rhs_ode, [0, Lstar], [0.0, 0.0, par[0], par[1]],
                   t_eval=[Lstar/2, Lstar], rtol=1e-12, atol=1e-15)
    return [so.y[2,1] - CM*Mstar, so.y[3,1]]
par = fsolve(shoot, [CM*Mstar, 0.0])
so = solve_ivp(rhs_ode, [0, Lstar], [0.0, 0.0, par[0], par[1]],
               t_eval=[Lstar/2, Lstar], rtol=1e-12, atol=1e-15)
print(f"ind. solve_ivp: d(L/2)={so.y[0,0]*1e6:.5f} um (vs {d_mid*1e6:.5f}); "
      f"d(L)={so.y[0,1]*1e6:.5f} um")

# 2) small sigma_b limit -> bending-only: delta_t = CM*M*L^2/2
Mtest = 0.1
print("sigma_b->0 check (M=0.1): delta_t exact vs CM*M*L^2/2:",
      delta_end(5e-3, Mtest)*1e6, "vs", 0.5*CM*Mtest*(5e-3)**2*1e6)

# 3) energy check: -dPi/dL == b*(J_far - int_0^{dt} t dd) = M^2 CM/2 - b*int t dd
def Pi(L, M, nq=1200):
    s, w = np.polynomial.legendre.leggauss(nq)
    s = 0.5*L*(s+1); w = 0.5*L*w
    d, _ = solution(L, M, s)
    t = sig*(1 - d/dc)
    K = np.clip(s[None, :] - s[:, None], 0, None) * w[None, :]
    I = K @ t                       # I(s) = int_s^L (xi-s) t(xi) dxi
    m = M - b*I                     # m1 = m2 = m (same for both arms)
    beam = np.sum(w*(m**2/(2*D1) + m**2/(2*D2)))
    spring = b*np.sum(w*(sig*(d - 0.5*d**2/dc)))
    th = np.sum(w*m*(1/D1 + 1/D2))  # theta1(L)+theta2(L)
    return beam + spring - M*th

def G_energy(L, M, hL=2e-6):
    return -(Pi(L+hL, M) - Pi(L-hL, M))/(2*hL)

for (Lt, Mt) in [(Lstar, Mstar), (0.7*Lstar, 0.85*Mstar), (1.3*Lstar, 1.1*Mstar)]:
    dt = delta_end(Lt, Mt)
    ge = G_energy(Lt, Mt)
    gj = Mt**2*CM/2 - b*bridging_work(dt)
    print(f"L={Lt*1e3:7.3f} mm, M={Mt:.5f}: -dPi/dL={ge:12.6f}, "
          f"M^2CM/2-b*int tdd={gj:12.6f}, rel.diff={abs(ge-gj)/abs(gj):.2e}")

# 4) initiation limit L->0: M -> sqrt(2 b Gamma0/C_M)
Minit = np.sqrt(2*b*G0/CM)
def M_of_L(L):
    d0 = delta_end(L, 0.0)          # affine in M: delta_t = d0 + slope*M
    d1 = delta_end(L, 1.0)
    slope = d1 - d0
    c, s2 = sig, -sig/(2*dc)
    A2 = CM/(2*b) - s2*slope**2
    A1 = -2*s2*slope*d0 - c*slope
    A0 = -(G0 + c*d0 + s2*d0**2)
    r = [x.real for x in np.roots([A2, A1, A0]) if abs(x.imag) < 1e-10 and x.real > 0]
    return min(r)
print(f"M_init(LEFM) = {Minit:.6f} N m; M(L) along R-curve:")
for L in [1e-6, 1e-4, 1e-3, 3e-3, Lstar]:
    MM = M_of_L(L)
    print(f"  L={L*1e3:8.4f} mm: M={MM:.6f} N m, delta_t={delta_end(L,MM)*1e6:9.4f} um, "
          f"J={MM**2*CM/(2*b):7.4f} J/m^2, R=Gamma0+int tdd="
          f"{G0+bridging_work(delta_end(L,MM)):.4f}")
