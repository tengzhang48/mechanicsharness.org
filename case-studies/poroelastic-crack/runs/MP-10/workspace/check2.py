import sympy as sp

nu, nu_u, alp, G, kap, sig = sp.symbols('nu nu_u alpha G kappa sigma', positive=True)
p = sp.symbols('p', positive=True)

print("=== symbolic re-derivation of V1_short_coeff and cross-checks ===\n")

# constitutive laws as given in the problem
e_yy, e_xx, s_xx, s_yy, s_kk = sp.symbols('e_yy e_xx sigma_xx sigma_yy sigma_kk')
eq_eyy = sp.Eq(2*G*e_yy, s_yy - nu*s_kk + alp*(1-2*nu)*p)   # trace of the e_ij law (indices 1,2)
zeta = (alp*(1-2*nu)*s_kk + alp**2*(1-2*nu)**2*p/(nu_u-nu))/(2*G)

# (1) undrained pore pressure at the crack face: sigma_kk = -2*sigma, sigma_yy = -1
p_star = sp.solve(sp.Eq(zeta.subs(s_kk, -2*sig), 0), p)[0]
print("(1) p* (undrained face pressure) =", sp.simplify(p_star))
B = 3*(nu_u-nu)/(alp*(1-2*nu))
print("    check p* = 2B*sigma/3 :", sp.simplify(p_star - 2*B*sig/3) == 0)

# (2) 1D skin: sigma_yy = -sigma fixed, e_xx held fixed; solve sigma_xx, then de_yy/dp
e_inf = sp.symbols('e_inf')
s_xx_sol = sp.solve(sp.Eq(2*G*e_inf, s_xx - nu*(s_xx+s_yy) + alp*(1-2*nu)*p), s_xx)[0]
s_kk_sub = s_xx_sol + s_yy
eyy_expr = sp.solve(eq_eyy.subs(s_kk, s_kk_sub), e_yy)[0]
coeff_p = sp.simplify(sp.diff(eyy_expr, p))
print("(2) d e_yy/d p (skin, any fixed e_xx) =", sp.simplify(coeff_p),
      "  == alpha(1-2nu)/(2G(1-nu)):", sp.simplify(coeff_p - alp*(1-2*nu)/(2*G*(1-nu))) == 0)

# (3) storage under fixed e_xx -> diffusivity equals the given c
#     (zeta law is in terms of sigma_kk and p only; sigma_kk(p) via sigma_xx above)
S_p = sp.simplify(sp.diff(zeta.subs(s_kk, s_kk_sub), p))
c_skin = sp.simplify(kap/S_p)
c_given = 2*kap*G*(1-nu)*(nu_u-nu)/(alp**2*(1-2*nu)**2*(1-nu_u))
print("(3) skin diffusivity == given c:", sp.simplify(c_skin - c_given) == 0)

# (4) face displacement: u_y(0) = -int_0^inf e_yy dy measured from far field =>
#     du = -coeff * int (p - p*) dy ; p - p* = -p* erfc(y/2sqrt(ct))
#     int_0^inf erfc(y/(2 sqrt(ct))) dy = 2 sqrt(ct)/sqrt(pi)
t = sp.symbols('t', positive=True)
du_face = sp.simplify(-coeff_p * (-p_star*2*sp.sqrt(t)/sp.sqrt(sp.pi)))   # per unit sqrt(c t)
print("(4) du_face/(sqrt(c t)) =", sp.simplify(du_face))

# (5) volume increment: crack length 2L, opening change 2*du per unit x;
#     A = G*dV/L^2 / sqrt(ct/L^2)   (c cancels: du already carries sqrt(c t))
L = sp.symbols('L', positive=True)
A_sym = sp.simplify(2*L*2*du_face*G/L**2 * L/sp.sqrt(t))
print("(5) A =", sp.simplify(A_sym),
      " == 8 sigma (nu_u-nu)/((1-nu)sqrt(pi)):",
      sp.simplify(A_sym - 8*sig*(nu_u-nu)/((1-nu)*sp.sqrt(sp.pi))) == 0)

# (6) numeric
vals = {nu: sp.Rational(2,10), nu_u: sp.Rational(4,10), alp: sp.Rational(89,100),
        G: 1, L: 1, sig: 1}
print("    numeric A =", sp.N(A_sym.subs(vals), 12))
print("    numeric p* =", sp.N(p_star.subs(vals), 12))
print("    V1_0 = pi(1-nu_u) =", sp.N((sp.pi*(1-nu_u)).subs(vals), 12))
print("    V1_inf = pi(1-nu) =", sp.N((sp.pi*(1-nu)).subs(vals), 12))
print("    V2_inf = 2*alpha*(1-2nu)/pi =", sp.N((2*alp*(1-2*nu)/sp.pi).subs(vals), 12))
