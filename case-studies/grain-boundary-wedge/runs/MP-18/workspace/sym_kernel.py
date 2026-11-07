import sympy as sp

z, y = sp.symbols('z y', positive=True)
S = 1/(z - y) - 1/(z + y) - 2*y*(z - y)/(z + y)**3
F = sp.integrate(S, y)
F = sp.simplify(F)
print("Antiderivative F(y) =", F)
# derivative check
print("dF/dy - S =", sp.simplify(sp.diff(F, y) - S))

# March kernel: from primitives P1+P2+P3 (c0=1,c1=0)
K_march = 1/(y - z) + 1/(y + z) - 4*z/(y + z)**2 + 4*z**2/(y + z)**3
Fmarch = sp.integrate(K_march, y)
print("\nmarch kernel antideriv =", sp.simplify(Fmarch))
print("matches log|y-z|+log|y+z| + 4z/(y+z)-2z^2/(y+z)^2 ?",
      sp.simplify(Fmarch - (sp.log(y - z) + sp.log(y + z) + 4*z/(y+z) - 2*z**2/(y+z)**2)))

# Is statement S(x,y) antisymmetric in some sense? S(x,y)+S(y,x)?
S2 = 1/(y - z) - 1/(y + z) - 2*z*(y - z)/(y + z)**3
print("\nS(z,y)+S(y,z) =", sp.simplify(S + S2))
print("S(z,y)-S(y,z) =", sp.simplify(S - S2))