import sympy as sp

# ======================================================================
# Far-behind wake state, Case A, careful.
# ======================================================================
# Strip 0<=x2<=H, crack faces at x2=0 and x2=H for x1<0 (traction free).
# Far behind (x1->-inf) everything independent of x1? No -- foundation
# relaxation makes x1-dependence; but far behind the fields should
# approach a bounded equilibrium of the relaxational ODE.
#
# The "far behind reached relaxed cracked wake": faces traction free and
# the film has slid/opened. Since the foundation pulls u->0, the relaxed
# state is the bounded solution of the model with faces free.
# Independently of x1, a bounded state with u=const satisfies equilibrium
# h*div(sigma) - k u = 0. For uniform u, div sigma=0 automatic if sigma
# uniform. Traction free on x2=0,H: sigma22=0, sigma12=0 there.
# With uniform fields, sigma22=0 everywhere and sigma12=0. Also eps22=0.
# Then u2 = const (set 0). eps11 = u1'? if uniform, eps=0 -> sigma=0.
# But prestress sigma = s0 + Dsig; traction free needs sigma22=0:
#   s0 + Ebar*nu*eps11=0 -> wait for uniform u1 = a constant then eps=0
#   -> sigma22=s0 !=0. So uniform u is NOT traction free. Need eps11 so
#   that s0 + Dsig_22 =0 and s0 + Dsig_11 =0 simultaneously if sigma22 and
#   sigma11 both zero? On faces only sigma22,sigma12 need vanish.
# Set plane stress, u1 = u1(x2)? Hmm.
#
# KEY: far behind in the WAKE the whole strip between two cracks relaxes
# uniaxially along x1: faces at x2=0,H free => sigma22=0 everywhere across
# thickness (since H small relative to variation across? no).
# Actually since the cracks run along x1 and are periodic in x2, far
# behind the solution is independent of x1: a 1D problem in x2.
# Let u = (u1(x2), u2(x2)).
# eps22 = u2', eps12 = u1'/2, eps11=0.
# traction free faces sigma22=0, sigma12=0 at x2=0,H.
x2=sp.symbols('x2',real=True)
nu=sp.symbols('nu',positive=True)
Ebar=sp.symbols('Ebar',positive=True)
s0=sp.symbols('sigma0',positive=True)
k=sp.symbols('k',positive=True)
h=sp.symbols('h',positive=True)
u1=sp.Function('u1');u2=sp.Function('u2')
eps11=0;eps22=sp.diff(u2(x2),x2);eps12=sp.Rational(1,2)*sp.diff(u1(x2),x2)
tr=eps11+eps22
sig11=Ebar*((1-nu)*eps11+nu*tr)
sig22=Ebar*((1-nu)*eps22+nu*tr)
sig12=Ebar*(1-nu)*eps12
# equilibrium: h sig_{ab,b}-k u_a=0
eq1=sp.simplify(h*sp.diff(sig12,x2)-k*u1(x2))
eq2=sp.simplify(h*sp.diff(sig22,x2)-k*u2(x2))
print("eq1:",eq1)
print("eq2:",eq2)
sol1=sp.dsolve(sp.Eq(eq1,0),u1(x2))
sol2=sp.dsolve(sp.Eq(eq2,0),u2(x2))
print("u1:",sol1)
print("u2:",sol2)