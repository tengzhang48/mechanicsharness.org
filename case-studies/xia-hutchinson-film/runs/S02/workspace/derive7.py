import sympy as sp

# ======================================================================
# Far-behind relaxed wake, 1D in x2, solve BVP exactly.
# u1(x2): u1'' = (2k/(Ebar h (1-nu))) u1 ; length L1 = sqrt(Ebar h (1-nu)/(2k))
# u2(x2): u2'' = (k/(Ebar h)) u2          ; length L2 = sqrt(Ebar h/k) = ell
# Faces x2=0,H traction free: sig22=0, sig12=0.
# ----------------------------------------------------------------------
x2,H,nu,Ebar,s0,k,h=sp.symbols('x2 H nu Ebar sigma0 k h',positive=True)
C1,C2,D1,D2=sp.symbols('C1 C2 D1 D2')
L1=sp.sqrt(Ebar*h*(1-nu)/(2*k))
L2=sp.sqrt(Ebar*h/k)      # = ell
# choose hyperbolic form to make BCs clean
u1=(C1*sp.cosh((x2-H/2)/L1)+C2*sp.sinh((x2-H/2)/L1))
u2=(D1*sp.cosh((x2-H/2)/L2)+D2*sp.sinh((x2-H/2)/L2))
eps11=0;eps22=sp.diff(u2,x2);eps12=sp.Rational(1,2)*sp.diff(u1,x2)
tr=eps22
sig11=Ebar*((1-nu)*eps11+nu*tr)
sig22=Ebar*((1-nu)*eps22+nu*tr)
sig12=Ebar*(1-nu)*eps12
# BC at x2=0 and H: sig22=0, sig12=0
eqs=[sig22.subs(x2,0),sig22.subs(x2,H),sig12.subs(x2,0),sig12.subs(x2,H)]
sol=sp.solve(eqs,[D1,D2,C1,C2],dict=True)
print(sol)
for s in sol:
    print("u1:",sp.simplify(u1.subs(s)))
    print("u2:",sp.simplify(u2.subs(s)))