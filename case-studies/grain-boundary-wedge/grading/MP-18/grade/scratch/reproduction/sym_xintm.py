import numpy as np
import sympy as sp

x, y = sp.symbols('x y', positive=True)
S = 1/(y - x) + 1/(x + y) + 2*y*(x - y)/(x + y)**3

# Antiderivatives Int x^m S dx for m=0..6, as expressions; then evaluate at 0,1.
exprs = []
for m in range(0, 8):
    Fm = sp.integrate(x ** m * S, x)
    exprs.append(sp.simplify(Fm))
    print("m=%d ok, d/dx check:" % m, sp.simplify(sp.diff(Fm, x) - x ** m * S))

# Build lambdified per-m expressions of the DEFINITE integral G_m(y) = Fm(1)-Fm(0)
import pickle
lams = []
for m in range(8):
    Gm = sp.simplify(exprs[m].subs(x, 1) - exprs[m].subs(x, 0))
    # sp.lambdify with numpy
    f = sp.lambdify(y, Gm, 'numpy')
    lams.append(f)
    print("G_%d(y) = %s" % (m, sp.simplify(Gm)))

# save lambdified funcs via a generated module
with open("xint_m.py", "w") as fh:
    fh.write("import numpy as np\n")
    for m in range(8):
        Gm = sp.simplify(exprs[m].subs(x, 1) - exprs[m].subs(x, 0))
        src = sp.printing.pycode(Gm)
        src = src.replace("math.", "np.")
        fh.write("def G%d(y):\n    return %s\n\n" % (m, src))
print("wrote xint_m.py")