import numpy as np
import hrr2 as H

# Linear-elastic limit n=1: p=1/2, k=3/2.  Williams Mode I angular functions.
H.n = 1.0
H.p = 0.5
H.k = 1.5

t = np.linspace(0.0, np.pi, 13)
sr = np.cos(t/2)*(1 - np.sin(t/2)*np.sin(3*t/2))
st = np.cos(t/2)*(1 + np.sin(t/2)*np.sin(3*t/2))
srt = np.cos(t/2)*np.sin(t/2)*np.cos(3*t/2)

# B = k(k-1) phi = 0.75 phi  => phi = st/0.75
# analytic: phi = st/0.75 = (4/3) cos(t/2)(1+sin(t/2)sin(3t/2))
import sympy as sp
tt = sp.symbols('t')

# Exact Airy eigenfunction for Mode I with p=1/2 (Williams lambda=1/2):
# phi = cos(3t/2) + (1/3) cos(t/2)  (the combination that is traction-free at t=pi)
phi_s = sp.cos(3*tt/2) + 3*sp.cos(tt/2)
d1_s = sp.diff(phi_s, tt)
d2_s = sp.diff(phi_s, tt, 2)
d3_s = sp.diff(phi_s, tt, 3)
A_s = sp.Rational(3, 2)*phi_s + d2_s          # sigma_rr
B_s = sp.Rational(3, 4)*phi_s                 # sigma_tt
C_s = -sp.Rational(1, 2)*d1_s                 # sigma_rt
fA = sp.lambdify(tt, sp.simplify(A_s), 'numpy')
fB = sp.lambdify(tt, sp.simplify(B_s), 'numpy')
fC = sp.lambdify(tt, sp.simplify(C_s), 'numpy')
print("Airy eigenfunction phi = cos(3t/2)+(1/3)cos(t/2):  check free surface & connector")
for x in [0.0, np.pi/2, np.pi]:
    print("  t=%.3f  A=%.5f  B=%.5f  C=%.5f" % (x, fA(x), fB(x), fC(x)))

# Validate the compatibility ODE at n=1 against this eigenfunction.
# Build state Y=[phi,phi',phi'',phi'''] from sympy derivatives, evaluate residual.
f1 = sp.lambdify(tt, d1_s, 'numpy')
f2 = sp.lambdify(tt, d2_s, 'numpy')
f3 = sp.lambdify(tt, d3_s, 'numpy')
f0 = sp.lambdify(tt, phi_s, 'numpy')
print("residual phi'''' - rhs(phi) at n=1:")
for x in [0.3, 1.0, 2.0, 2.8]:
    Y = [f0(x), f1(x), f2(x), f3(x)]
    r = H.rhs(x, Y)           # returns [d1,d2,d3,d3p]
    d4_analytic = sp.lambdify(tt, sp.diff(phi_s, tt, 4), 'numpy')(x)
    print("  t=%.2f  d4 analytic=%.5f  rhs d3p=%.5f  diff=%.2e" %
          (x, d4_analytic, r[3], d4_analytic - r[3]))