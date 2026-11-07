import numpy as np
import sympy as sp

# Generate antiderivatives F_0=S_phys, F_1=Int S dy, F_2=Int F_1 dy, ... up to order D.
# Then Int y^m S dy = y^m F_1 - m y^{m-1} F_2 + m(m-1) y^{m-2} F_3 - ... (repeated IBP).
# Specifically  Int y^m S dy = sum_{k=1}^{m+1} (-1)^{k-1} (m)_{k-1} y^{m-k+1} F_k.
x, y = sp.symbols('x y', positive=True)
S = 1/(y - x) + 1/(x + y) + 2*y*(x - y)/(x + y)**3

D = 16
Fs = [S]
for k in range(1, D + 1):
    Fs.append(sp.simplify(sp.integrate(Fs[-1], y)))

# verify recurrence for Int y^m S dy
def Int_ymS(m):
    import math
    tot = 0
    for k in range(1, m + 2):
        c = (-1) ** (k - 1) * (math.factorial(m) // math.factorial(m - k + 1))
        tot += c * y ** (m - k + 1) * Fs[k]
    return sp.simplify(tot)

# check against integration
mchk = 3
print("check m=3:", sp.simplify(sp.diff(Int_ymS(mchk), y) - y ** mchk * S))

# Generate lambdified DEFINITE integrals  K_m(x) = Int_0^1 y^m S dy
with open("ymint.py", "w") as fh:
    fh.write("import numpy as np\n")
    for m in range(0, D):
        Km = sp.simplify(Int_ymS(m).subs(y, 1) - Int_ymS(m).subs(y, 0))
        src = sp.printing.pycode(Km).replace("math.", "np.")
        fh.write("def K%d(x):\n    return %s\n\n" % (m, src))
    fh.write("KU = [")
    fh.write(", ".join("K%d" % m for m in range(0, D)))
    fh.write("]\n")
print("wrote ymint.py with K_0..K_%d" % (D - 1))
print("sample K_0:", sp.simplify(Int_ymS(0).subs(y, 1) - Int_ymS(0).subs(y, 0)))