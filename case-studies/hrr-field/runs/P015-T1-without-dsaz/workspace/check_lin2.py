import numpy as np
import hrr2 as H

# Find the true Mode I Airy eigenfunction for p=1/2: phi = c1 cos(3t/2)+c2 cos(t/2)
# solved so that free surface conditions A(pi)=0 and C(pi)=0 hold.
# A = k phi + phi'', k=3/2 ; C = -(k-1) phi' = -0.5 phi'.
import sympy as sp
t = sp.symbols('t')
c1, c2 = sp.symbols('c1 c2')
phi = c1*sp.cos(3*t/2) + c2*sp.cos(t/2)
A = sp.Rational(3, 2)*phi + sp.diff(phi, t, 2)
C = -sp.Rational(1, 2)*sp.diff(phi, t)
A_pi = sp.simplify(A.subs(t, sp.pi))
C_pi = sp.simplify(C.subs(t, sp.pi))
print("A(pi) =", A_pi)
print("C(pi) =", C_pi)
sol = sp.solve([A_pi, C_pi], [c1, c2], dict=True)
print("solve:", sol)

# pick c1=1 -> c2?
expr = sp.solve(A_pi, c1)
print("c1 in terms of c2:", expr)