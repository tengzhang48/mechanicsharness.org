"""
Clean derivation + symbolic consistency check of the HRR angular equation.

Airy:  Phi = r^k phi(theta), k = 2 - p, p = 1/(n+1).
Stress (normalized, factor r^{-p}):
  A = sigtilde_rr = k phi + phi''
  B = sigtilde_tt = k(k-1) phi
  C = sigtilde_rt = -(k-1) phi'
Equilibrium automatic.

Plane strain: sig_zz=(sig_rr+sig_tt)/2 => s_zz=0.
  s_rr=(A-B)/2, s_tt=(B-A)/2, s_rt=C.
  se^2 = (3/4)(A-B)^2 + 3 C^2.

Constitutive: eps = (3/2) alpha eps0 (se/s0)^{n-1} s/s0.
Displacement: u_r = C0 r^{-p} Ur, u_t = C0 r^{-p} Ut, C0 = alpha eps0 K^n.
  eps_rr = C0 r^{-p-1}(-p Ur)
  eps_tt = C0 r^{-p-1}(Ur + Ut')
  2 eps_rt = C0 r^{-p-1}(Ur' - p Ut - Ut)
Constitutive (q=se^{n-1}):
  -p Ur = q (A-B)/2
  Ur + Ut' = q (B-A)/2
  Ur' - p Ut - Ut = 2 q C
"""
import sympy as sp

t = sp.symbols('t')
p = sp.symbols('p', positive=True)
nn = sp.symbols('n', positive=True)
k = 2 - p

phi = sp.Function('phi')(t)
A = k*phi + sp.diff(phi, t, 2)
B = k*(k-1)*phi
C = -(k-1)*sp.diff(phi, t)

se2 = sp.Rational(3, 4)*(A-B)**2 + 3*C**2
se = sp.sqrt(se2)
q = se**(nn-1)

Ur = -q*(A-B)/(2*p)          # from eq1
Utp_relation = q*(B-A)/2     # = Ur + Ut'  => Ut' = this - Ur
Utp = Utp_relation - Ur
Urp = sp.diff(Ur, t)
Ut = (Urp - 2*q*C)/(p+1)     # from eq3
eq = sp.expand(sp.diff(Ut, t) - Utp)

d4 = sp.diff(phi, t, 4)
sol = sp.solve(sp.Eq(eq, 0), d4)
print("Got closed form for phi'''' =", sol)