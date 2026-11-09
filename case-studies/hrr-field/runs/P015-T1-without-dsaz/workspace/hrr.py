"""
HRR plane-strain Mode I crack-tip field for power-law hardening material.

Constitutive (total deformation):
  eps_ij = (3/2) alpha eps0 (sigma_e/sigma0)^(n-1) s_ij / sigma0
  s_ij = sigma_ij - (1/3) sigma_kk delta_ij,  sigma_e = sqrt(3/2 s_ij s_ij)

Self-similar leading field:
  sigma_ij(r,theta) = sigma0 * K * r^{-p} * sigtilde_ij(theta),  p = 1/(n+1)

Angular stress functions (normalized by sigma0*K):
  sigma_rr = r^{-p} a(theta)
  sigma_tt = r^{-p} b(theta)
  sigma_rt = r^{-p} c(theta)

Equilibrium gives:
  c' = a (p-1) + b
  b' = c (p-2)

Plane strain: eps_zz=0 => s_zz=0 => sigma_zz=(param).
  s_rr = (a-b)/2, s_tt = (b-a)/2, s_rt = c
  se^2 = (3/4)(a-b)^2 + 3 c^2

Displacement: u_r = r^{-p} A(theta), u_t = r^{-p} B(theta)
  eps_rr = -p r^{-p-1} A
  eps_tt = r^{-p-1} (A + B')
  2 eps_rt = r^{-p-1} (A' - pB - B)
Constitutive matching (common factor Lambda r^{-p-1}):
  -p A = se^{n-1} s_rr
  A + B' = se^{n-1} s_tt
  A' - pB - B = 2 se^{n-1} s_rt

=> A = -q u/(2p),  u = a-b,  q = se^{n-1}
   B' = q u (1-p)/(2p)
   A' = 2 q c + (p+1) B
Compatibility: differentiate A, and enforce A'=... gives an ODE for a (2nd order).

We integrate as a first-order system with shooting.
"""
import numpy as np
from scipy.integrate import solve_ivp, solve_bvp
from scipy.optimize import brentq

n = 6.5
alpha = 0.02
eps0 = 400e6/200e9
p = 1.0/(n+1.0)
print("n =", n, " alpha =", alpha, " eps0 =", eps0)
print("p = 1/(n+1) =", p)


# ---------------------------------------------------------------------------
# Anglar ODE system (HRR). State y = [a, b, c, ap] where ap = a'(theta).
#   b'  = c (p-2)
#   c'  = a (p-1) + b
#   a'  = ap
#   ap' = from compatibility relation below
# ---------------------------------------------------------------------------
def se_of(a, b, c):
    return np.sqrt(0.75*(a-b)**2 + 3.0*c**2)


def rhs(theta, y):
    a, b, c, ap = y
    # derivatives of b, c
    bp = c*(p-2.0)
    cp = a*(p-1.0) + b
    u = a - b
    se = se_of(a, b, c)
    q = se**(n-1.0)
    # second derivative of b
    bpp = cp*(p-2.0)
    # u' and se'
    up = ap - bp
    sep = (1.5*u*up + 6.0*c*cp)/(2.0*se)
    qp = (n-1.0)*se**(n-2.0)*sep
    # A = -q u/(2p);  A' = 2 q c + (p+1) B  => B = (A' - 2 q c)/(p+1)
    A = -q*u/(2.0*p)
    Ap = -(1.0/(2.0*p))*(qp*u + q*up)
    B = (Ap - 2.0*q*c)/(p+1.0)
    # B' = (A'' - 2(q' c + q c'))/(p+1)  must equal q u (1-p)/(2p)
    # Solve for A'' then get a'':
    #   A'' = 2(q' c + q c') + (p+1) q u (1-p)/(2p)
    App = 2.0*(qp*c + q*cp) + (p+1.0)*q*u*(1.0-p)/(2.0*p)
    # A'' = -(1/(2p))[ q'' u + 2 q' u' + q (a''-b'') ]
    # => q'' u + 2 q' u' + q (app - bpp) = -2p App
    # need q'':  sep = (1.5 u up + 6 c cp)/(2 se); sep2 = derivative
    # se' = (1.5 u up + 6 c cp)/(2 se).  Compute se'' numerically-consistently:
    # Use identity se se'' + sep^2 = (1/2) d2(se^2)/dth2 with se^2=(3/4)u^2+3c^2
    # (se^2)'' = (3/2)(u'^2+u u'') + 6(c'^2 + c c'') ; u''=app-bpp, c''=cp'
    cpp = ap*(p-1.0) + bp   # c'' = a'(p-1)+b'
    # Solve implicitly for L = a''. Collect L-dependent terms:
    #   se2pp = (3/2) u L + R0,   R0 = (3/2)u'^2 + 6 c'^2 + 6 c c'' - (3/2) u b''
    R0 = 1.5*up*up + 6.0*cp*cp + 6.0*c*cpp - 1.5*u*bpp
    #   se'' = (3/4)(u/se) L + (R0 - 2 se'^2)/(2 se)
    #   q'' = (n-1)se^{n-2} se'' + (n-1)(n-2) se^{n-3} se'^2
    #       = alpha_q L + beta_q
    alpha_q = (n-1.0)*se**(n-2.0)*0.75*u/se          # se^{n-2}*(3/4)(u/se) = (3/4) se^{n-3} u
    beta_q = (n-1.0)*se**(n-2.0)*(R0 - 2.0*sep*sep)/(2.0*se) \
             + (n-1.0)*(n-2.0)*se**(n-3.0)*sep*sep
    # Equation: (alpha_q L + beta_q) u + 2 qp up + q (L - bpp) = -2 p App
    num = -2.0*p*App - beta_q*u - 2.0*qp*up + q*bpp
    den = alpha_q*u + q
    app = num/den
    return [ap, bp, cp, app]


# Boundary conditions for solve_bvp.
# State y = [a, b, c, ap];  theta in [0, pi].
# theta=0 axis (Mode I symmetry):  c(0)=0, ap(0)=0
# theta=pi free surface:          a(pi)=0, c(pi)=0
# Scaling: b(0)=1 fixes normalization.
def bc(ya, yb):
    return np.array([ya[2], ya[3], yb[0], yb[2]])


def solve_hrr():
    th = np.linspace(0.0, np.pi, 400)
    # initial guess: smooth functions satisfying end conditions roughly
    a0 = a0_guess
    y0 = np.zeros((4, th.size))
    y0[0] = a0*np.ones_like(th)             # a
    y0[1] = np.ones_like(th)                # b
    y0[2] = 0.0                             # c
    y0[3] = 0.0                             # ap
    sol = solve_bvp(rhs, bc, th, y0, tol=1e-10, max_nodes=20000, verbose=0)
    return sol


if __name__ == "__main__":
    a0_guess = 0.3
    sol = solve_hrr()
    print("solve_bvp success:", sol.success, " message:", sol.message)
    print("nodes:", sol.x.size)