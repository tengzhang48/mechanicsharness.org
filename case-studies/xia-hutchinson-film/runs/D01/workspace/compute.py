import sympy as sp

# Symbols
x, H, ell, s0, Ebar = sp.symbols('x H ell sigma0 Ebar', positive=True)
nu, h, k = sp.symbols('nu h k', positive=True)

# ---------------------------------------------------------------------------
# 1. Governing 2D reduction (symbolic check of the PDE)
# ---------------------------------------------------------------------------
u1, u2 = sp.Function('u1')(x), sp.Function('u2')(x)
E11, E22 = sp.diff(u1, x)*0, sp.diff(u2, x)   # only x2-dependence => e11=0
# in-plane strains (x1-invariant): e11 = 0, e22 = u2', e12=(u1')/2
trE = E22
# equilibrium component equations (coefficient of Ebar after division):
# alpha=1: (1-nu)/2 u1'' - k/(h Ebar) u1 = 0
# alpha=2: u2''        - k/(h Ebar) u2 = 0
eq1 = sp.Eq((1-nu)/2*sp.diff(u1, x, 2) - k/(h*Ebar)*u1, 0)
eq2 = sp.Eq(sp.diff(u2, x, 2) - k/(h*Ebar)*u2, 0)

# wake length: 1/ell^2 = k/(h Ebar)
ell_expr = sp.sqrt(h*Ebar/k)
print("ell =", ell_expr)

# Substitute k/(h Ebar) = 1/ell^2 for the wake ODEs
eq2w = eq2.subs(k/(h*Ebar), 1/ell**2)
eq1w = eq1.subs(k/(h*Ebar), 1/ell**2)

# Solve wake ODEs
sol2 = sp.dsolve(eq2, u2)
sol1 = sp.dsolve(eq1, u1)
print("u2 general:", sol2.rhs)
print("u1 general:", sol1.rhs)

# ---------------------------------------------------------------------------
# 2. Case A far-behind wake, strip width H, edges at x2=0,H traction free
# ---------------------------------------------------------------------------
A_, B_ = sp.symbols('A B')
u2w = A_*sp.cosh((x - H/2)/ell) + B_*sp.sinh((x - H/2)/ell)
# BC: sigma22 = sigma0 + Ebar*u2' = 0 at x=0 and x=H
bc0 = sp.Eq(s0 + Ebar*sp.diff(u2w, x).subs(x, 0), 0)
bcH = sp.Eq(s0 + Ebar*sp.diff(u2w, x).subs(x, H), 0)
solAB = sp.solve([bc0, bcH], [A_, B_], dict=True)[0]
print("A,B:", solAB)
u2A = sp.simplify(u2w.subs(solAB))
print("u2 wake (case A):", u2A)
print("u2'(0):", sp.simplify(sp.diff(u2A, x).subs(x, 0)))
print("u2'(H):", sp.simplify(sp.diff(u2A, x).subs(x, H)))
print("sigma22(x):", sp.simplify(s0 + Ebar*sp.diff(u2A, x)))

# energy density per unit x1-length (case A), width H
u2p = sp.diff(u2A, x)
w_dens = h*s0*u2p + sp.Rational(1,2)*h*Ebar*(u2p**2) + sp.Rational(1,2)*k*u2A**2
WA = sp.simplify(sp.integrate(w_dens, (x, 0, H)))
WA = sp.simplify(WA.subs(k, h*Ebar/ell**2))
print("W_behind case A =", sp.simplify(WA))
WA_expected = -h*s0**2*ell/Ebar*sp.tanh(H/(2*ell))
print("W_A - expected:", sp.simplify(WA - WA_expected))

# ---------------------------------------------------------------------------
# 3. Case B: period 2H, old cracks at x2=0,2H ; new crack at x2=H
# ---------------------------------------------------------------------------
# far ahead: strip width 2H
WA_head = WA_expected.subs(H, 2*H)          # per period 2H (width 2H)
# far behind: two strips width H
WB_behind = 2*WA_expected                    # two strips of width H
rel_B = sp.simplify(WA_head - WB_behind)     # release per period per unit advance
GB = sp.simplify(rel_B/h)
GA = sp.simplify((-WA_expected)/h)
print("GA =", GA)
print("GB =", GB)
print("gamma_A =", sp.simplify(GA*Ebar/(s0**2*ell)))
print("gamma_B =", sp.simplify(GB*Ebar/(s0**2*ell)))

# ---------------------------------------------------------------------------
# 4. Numeric sampling
# ---------------------------------------------------------------------------
import json
nu_v = sp.Rational(3,10)
vals = []
for eH in [sp.Rational(1,2), 1, sp.Rational(3,2), 2]:
    beta = sp.Rational(1,2)/eH          # beta = H/(2ell) = (1/(2))*(H/ell); H/ell=1/eH
    beta = (1/(2*eH))
    gA = sp.tanh(beta)
    gB = 2*sp.tanh(beta) - sp.tanh(2*beta)
    kIA = sp.sqrt((1-nu_v**2)*gA)
    kIB = sp.sqrt((1-nu_v**2)*gB)
    row = {
        "ell_over_H": float(eH),
        "gamma_A": float(sp.N(gA, 20)),
        "gamma_B": float(sp.N(gB, 20)),
        "kI_A": float(sp.N(kIA, 20)),
        "kI_B": float(sp.N(kIB, 20)),
    }
    vals.append(row)
    print(row)

with open("results_analytic.json", "w") as f:
    json.dump({"note": "analytic", "rows": vals}, f, indent=2)

# ---------------------------------------------------------------------------
# 5. Limits
# ---------------------------------------------------------------------------
beta = sp.symbols('beta', positive=True)
print("gamma_A limit beta->0 :", sp.limit(sp.tanh(beta), beta, 0))
print("gamma_A limit beta->oo:", sp.limit(sp.tanh(beta), beta, sp.oo))
print("gamma_B limit beta->0 :", sp.series(2*sp.tanh(beta)-sp.tanh(2*beta), beta, 0, 4))
print("gamma_B limit beta->oo:", sp.limit(2*sp.tanh(beta)-sp.tanh(2*beta), beta, sp.oo))
