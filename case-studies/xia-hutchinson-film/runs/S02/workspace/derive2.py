import sympy as sp

p, t = sp.symbols('p t', positive=True)
integrand = t*(1-t)/sp.sqrt((t+p)**2 - 1)
I = sp.integrate(integrand, (t, 0, 1))
print("I =", sp.simplify(I))
print()
F = p - 1 - I
print("F =", sp.simplify(F))
print()
for pv in [0.5,1.0,2.0,5.0,10.0]:
    print(pv, float(F.subs(p, pv)))