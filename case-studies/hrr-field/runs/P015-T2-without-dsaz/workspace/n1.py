import sympy as sp
th, s = sp.symbols('theta s')
U = sp.Function('U')(th); V = sp.Function('V')(th); P = sp.Function('P')(th)
# n=1 incompressible: eps_ij = s_ij (B=1). shapes:
d_rr = (1-s)*U
d_tt = U + sp.diff(V, th)
d_rt = sp.Rational(1,2)*(sp.diff(U, th) - s*V)
# stresses (deviatoric) shape + pressure P
S_rr = d_rr + P; S_tt = d_tt + P; S_rt = d_rt
eqr = sp.simplify(-s*S_rr + (S_rr - S_tt) + sp.diff(S_rt, th))       # r-eq
eqt = sp.simplify(-s*S_rt + 2*S_rt + sp.diff(S_tt, th))              # theta-eq
inc = sp.simplify((1-s)*U + U + sp.diff(V, th))                      # incompressibility
print("eqr =", sp.expand(eqr))
print("eqt =", sp.expand(eqt))
print("inc =", sp.expand(inc))
# eliminate V: from inc, V' = -s U ; V'' = -s U'
ceq = sp.simplify(eqt - sp.diff(inc, th))  # combine to remove V''
print("eqt - d/dth(inc) =", sp.expand(ceq))
# so ceq involves U'' and P' ; eqr involves U' and V' and P'
# Use: from inc V'=-sU. Substitute into eqr.
eqr2 = sp.expand(eqr.subs(sp.diff(V,th), -s*U))
# eqt2: substitute V' and V''.
eqt2 = sp.expand(eqt.subs({sp.diff(V,th): -s*U, sp.diff(V,(th,2)): -s*sp.diff(U,th)}))
print("eqr2 =", eqr2)
print("eqt2 =", eqt2)