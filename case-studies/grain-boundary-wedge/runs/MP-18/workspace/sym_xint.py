import numpy as np
import sympy as sp

# Exact integrals  I_m(y) = Int_0^1 x^m S_phys(x,y) dx,  m=0..Nc.
# S_phys(x,y) = 1/(y-x) + 1/(x+y) + 2y(x-y)/(x+y)^3 .
# We integrate in x with y a parameter.  The PV is at x=y (inside (0,1)); the
# antiderivative handles log|y-x| whose jump gives the PV value.

x, y = sp.symbols('x y', positive=True)
S = 1/(y - x) + 1/(x + y) + 2*y*(x - y)/(x + y)**3
print("S check integrate w.r.t x once:")
Fx = sp.integrate(S, x)
print("Int S dx =", sp.simplify(Fx))
print("check d/dx:", sp.simplify(sp.diff(Fx, x) - S))