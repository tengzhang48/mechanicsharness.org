import sympy as sp

# ======================================================================
# RIGOROUS DERIVATION, CASE A
# ======================================================================
# Plane stress, Mode I. Let w = u_2. Cracks on x2 = mH, x1<0.
# Take one strip 0<=x2<=H. By symmetry about x2=0 and x2=H/2? No:
# the crack line is at x2=0 (and x2=H). Periodic array.
#
# Governing eq (mode I, u1 = 0 everywhere by symmetry? NOT exactly, but
# for a symmetric (mode I) problem with a periodic array we can take
# u1 = 0 and w=w(x1,x2)).
#
# Actually for a straight crack array under remote equibiaxial stress the
# mode I field has u1 even, puts KII=0. Keep u1 too, but Mode I dominant.
#
# Governing PDE (both components), separable in x1:  e^{s x1}:
#   plane stress operator L(u) - (k/h) u = 0
# Fourier in x2 with wavenumber q = 2 pi m / H -> operator:
#   Ebar * [ ((1-nu) s^2 - q^2(1-nu)/2 ...)] complex.
# Let's build the full 2x2 symbol matrix.

s, q = sp.symbols('s q', real=True)
nu = sp.symbols('nu')
Ebar = sp.symbols('Ebar', positive=True)
kh = sp.symbols('kh', positive=True)   # k/h

# Plane stress constitutive: sigma_ab = Ebar[(1-nu)eps_ab+nu eps_gg d_ab]
# Equilibrium: h sigma_{ab,b} - k u_a = 0 -> sigma_{ab,b} = (k/h) u_a
# Use u = (A,B) e^{s x1 + i q x2}
# eps_11 = s A ; eps_22 = i q B ; eps_12 = (s B + i q A)/2
# sigma_11 = Ebar[(1-nu) s A + nu (s A + i q B)]
# sigma_22 = Ebar[(1-nu) i q B + nu (s A + i q B)]
# sigma_12 = Ebar (1-nu) * (s B + i q A)/2   [since nu*eps_gg term off-diagonal=0? 
#   sigma_12 = Ebar[(1-nu) eps_12 + nu eps_gg delta_12]; delta_12=0]
# eq1: s sigma_11 + i q sigma_12 = kh A
# eq2: s sigma_12 + i q sigma_22 = kh B
A,B = sp.symbols('A B')
eps11 = s*A; eps22 = sp.I*q*B; eps12 = (s*B + sp.I*q*A)/2
tr = eps11+eps22
sig11 = Ebar*((1-nu)*eps11 + nu*tr)
sig22 = Ebar*((1-nu)*eps22 + nu*tr)
sig12 = Ebar*(1-nu)*eps12
eq1 = sp.expand(s*sig11 + sp.I*q*sig12 - kh*A)
eq2 = sp.expand(s*sig12 + sp.I*q*sig22 - kh*B)
M = sp.Matrix([[sp.expand(eq1/A if A!=0 else eq1), 0],[0, sp.expand(eq2/B)]])
# Build matrix properly by coefficients
e1A = sp.expand(eq1).coeff(A); e1B = sp.expand(eq1).coeff(B)
e2A = sp.expand(eq2).coeff(A); e2B = sp.expand(eq2).coeff(B)
Mat = sp.Matrix([[e1A, e1B],[e2A, e2B]])
print("Symbol matrix M(s,q):")
sp.pprint(sp.simplify(Mat))
det = sp.simplify(Mat.det())
print("\ndet =", sp.factor(det))