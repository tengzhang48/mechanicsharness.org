import sympy as sp

# Setup: plane stress, isotropic. ODE from h sigma_ab,b - k u_a = 0.
# Let Ebar = E/(1-nu^2). Plane stress isotropic:
# sigma = Ebar[(1-nu) eps + nu tr(eps) I].
x1, x2 = sp.symbols('x1 x2', real=True)
nu = sp.symbols('nu', positive=True)
E = sp.symbols('E', positive=True)
k = sp.symbols('k', positive=True)
h = sp.symbols('h', positive=True)
s0 = sp.symbols('sigma0', positive=True)
Ebar = E/(1-nu**2)

# Try a rich ansatz for far-behind. We'll check div form for candidate fields.
# Candidate: u1 = U1(x1), u2 = 0
U1 = sp.Function('U1')
# eps11 = U1', eps22=0, eps12=0
eps11 = sp.diff(U1(x1), x1)
eps22 = 0
eps12 = 0
tr = eps11
s11 = Ebar*((1-nu)*eps11 + nu*tr)
s22 = Ebar*((0) + nu*tr)
s12 = 0
# equilibrium x1: h ds11/dx1 - k U1 = 0 ; x2: h ds12/dx1 - k*0 = 0 -> 0
eq1 = sp.simplify(h*sp.diff(s11, x1) - k*U1(x1))
print("eq1 =", eq1)
# h*Ebar*(1-nu)*U1'' - k U1 = 0
sol = sp.dsolve(sp.Eq(h*Ebar*(1-nu)*sp.diff(U1(x1), x1, 2) - k*U1(x1), 0))
print("sol U1:", sol)