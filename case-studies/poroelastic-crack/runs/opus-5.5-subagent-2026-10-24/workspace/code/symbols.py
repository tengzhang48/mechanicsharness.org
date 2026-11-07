"""Derive the half-plane (y>0) Fourier symbols of plane-strain Biot poroelasticity
in the Laplace domain.  Fields ~ exp(i k x), k>0 here.
Solution = undrained elastic part (Papkovich, nu_u) + irrotational diffusive part.
BCs at y=0: sigma_xy=0, u_y=U, q_y=-kappa dp/dy = Q.  Output sigma_yy(0), p(0)."""
import sympy as sp

x, y = sp.symbols('x y', real=True)
k, s, G, c = sp.symbols('k s G c', positive=True)
nu, nuu, al = sp.symbols('nu nu_u alpha', positive=True)
A, B, D, U, Q = sp.symbols('A B D U Q')

eta = al*(1-2*nu)/(2*(1-nu))
M = 2*G*(nuu-nu)/(al**2*(1-2*nu)*(1-2*nuu))
lam = 2*G*nu/(1-2*nu)
lamu = 2*G*nuu/(1-2*nuu)
S = al**2*(1-2*nu)**2*(1-nuu)/(2*G*(nuu-nu)*(1-nu))
kap = c*S
m = sp.sqrt(k**2 + s/c)

E = sp.exp(sp.I*k*x - k*y)
psi0 = A*E; psiy = B*E
ux_u = (sp.diff(psi0 + y*psiy, x))/(2*G)
uy_u = (sp.diff(psi0 + y*psiy, y) - 4*(1-nuu)*psiy)/(2*G)
phi = D*sp.exp(sp.I*k*x - m*y)
ux_d = sp.diff(phi, x); uy_d = sp.diff(phi, y)

def strains(ux, uy):
    exx = sp.diff(ux, x); eyy = sp.diff(uy, y)
    exy = (sp.diff(ux, y) + sp.diff(uy, x))/2
    return exx, eyy, exy

exx_u, eyy_u, exy_u = strains(ux_u, uy_u)
e_u = exx_u + eyy_u
sxx_u = 2*G*exx_u + lamu*e_u; syy_u = 2*G*eyy_u + lamu*e_u; sxy_u = 2*G*exy_u
p_u = -al*M*e_u
exx_d, eyy_d, exy_d = strains(ux_d, uy_d)
e_d = exx_d + eyy_d
p_d = (G/eta)*(s/c)*phi
sxx_d = 2*G*exx_d + lam*e_d - al*p_d
syy_d = 2*G*eyy_d + lam*e_d - al*p_d
sxy_d = 2*G*exy_d

# checks: undrained Navier, diffusive: equilibrium and mass balance
def div_s(sxx, syy, sxy):
    return [sp.simplify(sp.diff(sxx, x) + sp.diff(sxy, y)), sp.simplify(sp.diff(sxy, x) + sp.diff(syy, y))]
print('equil undrained:', div_s(sxx_u, syy_u, sxy_u))
print('equil diffusive:', div_s(sxx_d, syy_d, sxy_d))
zeta_u = al*e_u + p_u/M
print('zeta undrained:', sp.simplify(zeta_u), ' lap p_u:', sp.simplify(sp.diff(p_u, x, 2) + sp.diff(p_u, y, 2)))
zeta_d = al*e_d + p_d/M
mass_d = sp.simplify(s*zeta_d - kap*(sp.diff(p_d, x, 2) + sp.diff(p_d, y, 2)))
print('mass diffusive:', mass_d)

sxx = sxx_u + sxx_d; syy = syy_u + syy_d; sxy = sxy_u + sxy_d
uy = uy_u + uy_d; p = p_u + p_d
q = -kap*sp.diff(p, y)
sub0 = {y: 0}
eqs = [sp.Eq((sxy/sp.exp(sp.I*k*x)).subs(sub0), 0),
       sp.Eq((uy/sp.exp(sp.I*k*x)).subs(sub0), U),
       sp.Eq((q/sp.exp(sp.I*k*x)).subs(sub0), Q)]
sol = sp.solve(eqs, [A, B, D], dict=True)[0]
syy0 = sp.simplify((syy/sp.exp(sp.I*k*x)).subs(sub0).subs(sol))
p0 = sp.simplify((p/sp.exp(sp.I*k*x)).subs(sub0).subs(sol))
A11 = sp.simplify(sp.diff(syy0, U)); A12 = sp.simplify(sp.diff(syy0, Q))
A21 = sp.simplify(sp.diff(p0, U)); A22 = sp.simplify(sp.diff(p0, Q))
for nm, ex in [('A11', A11), ('A12', A12), ('A21', A21), ('A22', A22)]:
    print(nm, '=', sp.factor(sp.simplify(ex)))
import pickle
with open('<opus-workspace>/code/symbols.pkl', 'wb') as f:
    pickle.dump({'A11': A11, 'A12': A12, 'A21': A21, 'A22': A22}, f)
