import sympy as sp

z, y = sp.symbols('z y', positive=True)
S = 1/(z - y) - 1/(z + y) - 2*y*(z - y)/(z + y)**3

F1 = sp.simplify(sp.integrate(S, y))           # Int S dy
F2 = sp.simplify(sp.integrate(F1, y))           # Int F1 dy  (for Int y S dy)

print("F1 =", F1)
print()
print("F2 =", sp.simplify(F2))
print()
# Int (alpha + beta y) S dy = alpha*F1 + beta*(y*F1 - F2)
# verify
alpha, beta = sp.symbols('alpha beta')
expr = alpha*F1 + beta*(y*F1 - F2)
print("d/dy [alpha F1 + beta(y F1 - F2)] - (alpha+beta y)S =",
      sp.simplify(sp.diff(expr, y) - (alpha + beta*y)*S))

# print F2 in a simple form for numpy translation
print("\nF2 simplified:", sp.simplify(sp.expand(F2)))
print("\nF1 expanded:", sp.expand(F1))