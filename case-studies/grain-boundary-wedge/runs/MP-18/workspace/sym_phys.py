import numpy as np
import sympy as sp

# ============ Definitive kernel (problem statement, fixed sign) ============
# Statement:  s = 1 - (1/pi) PV Int_0^1 S(x,y) w(y) dy ,
#             S(x,y) = 1/(x-y) - 1/(x+y) - 2y(x-y)/(x+y)^3 .
# Physics requires decay; the singular term must act like +1/(y-x) on w
# (positive opening reduces traction).  We verified numerically that the
# decaying branch corresponds to  -S_statement  (overall sign flip), i.e.
#   S_phys(x,y) = 1/(y-x) + 1/(x+y) + 2y(x-y)/(x+y)^3 .
# We build  M_ij = (1/pi) Int_0^1 S_phys(x_i,y) phi_j(y) dy  EXACTLY via the
# sympy-derived antiderivative.
#
# Antiderivative of S_phys w.r.t. y  (sympy):
z, y = sp.symbols('z y', positive=True)
S_phys = 1/(y - z) + 1/(z + y) + 2*y*(z - y)/(z + y)**3
F1 = sp.simplify(sp.integrate(S_phys, y))
F2 = sp.simplify(sp.integrate(F1, y))
print("F1 =", F1)
print("F2 =", F2)
print("check:", sp.simplify(sp.diff(F1, y) - S_phys))
print("check2:", sp.simplify(sp.diff(F2, y) - F1))