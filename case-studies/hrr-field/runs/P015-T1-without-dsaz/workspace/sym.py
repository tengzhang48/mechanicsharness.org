import sympy as sp
th=sp.symbols('theta')
p,k=sp.symbols('p k', positive=True)
phi=sp.Function('phi')(th)
d1=sp.diff(phi,th); d2=sp.diff(phi,th,2); d3=sp.diff(phi,th,3); d4=sp.diff(phi,th,4)
A=k*phi+d2; B=k*(k-1)*phi; C=-(k-1)*d1
def dd(e): return sp.diff(e,th)
Cp=dd(C)
Ur=sp.Rational(3,4)*(A-B)/p          # n=1 => q=1
Urpp=dd(dd(Ur))
# eq: Ur'' = 3 C' - (p-1)[(3/4)(B-A) - Ur]
EQ=sp.Eq(Urpp, 3*Cp - (p-1)*(sp.Rational(3,4)*(B-A)-Ur))
sol=sp.solve(EQ, d4)
print("num sol:", len(sol))
print("d4 (n=1) =", sp.factor(sp.simplify(sol[0])))