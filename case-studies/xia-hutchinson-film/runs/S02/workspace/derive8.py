import sympy as sp

# ======================================================================
# Far-behind relaxed wake (x1-independent), strip [0,H], faces free.
# ======================================================================
x2,H,nu,Ebar,s0,k=sp.symbols('x2 H nu Ebar sigma0 k',positive=True)
a0,a1,b0,b1=sp.symbols('a0 a1 b0 b1')
# x1-independent: u1=u1(x2), u2=u2(x2)
# eps11=0, eps22=u2', eps12=u1'/2
u1=a0+a1*(x2-H/2)
u2=b0+b1*(x2-H/2)
eps11=0; eps22=sp.diff(u2,x2); eps12=sp.Rational(1,2)*sp.diff(u1,x2)
tr=eps22
sig22=s0+Ebar*((1-nu)*eps22+nu*tr)
sig12=Ebar*(1-nu)*eps12
# equilibrium: h*(d sig12/dx2) - k u1=0 ; h*(d sig22/dx2)-k u2=0
h=sp.symbols('h',positive=True)
eq1=h*sp.diff(sig12,x2)-k*u1
eq2=h*sp.diff(sig22,x2)-k*u2
# sig12, sig22 constant in x2 -> equilibrium gives -k u1=0,-k u2=0 => u=0.
print("eq1",sp.expand(eq1))
print("eq2",sp.expand(eq2))
# So with the foundation, the ONLY x1-independent bounded state is u=0,
# i.e. sigma=s0 everywhere (reference). So far behind, u->0 and sigma->s0 I.
# That means faces are NOT traction free far behind unless sig22=0 there.
# Contradiction => far-behind state must still vary in x1 (relax over ell).
print()
print("Conclusion: no x1-independent relaxed state. The wake decays over ell;")
print("the true far-behind limit of the DISPLACEMENT INCREMENT that matters is")
print("the bounded easterly decay; the crack faces are free LOCALLY, and the")
print("stress sig11 ahead->s0, behind->?")