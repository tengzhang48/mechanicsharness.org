import numpy as np
import sympy as sp

# Displacement formulation. u_i = r^{1-sn} U_i(theta). Incompressible power-law.
th, s, n = sp.symbols('theta s n')
U = sp.Function('U')(th); V = sp.Function('V')(th); P = sp.Function('P')(th)
r = sp.symbols('r', positive=True)

ur = r**(1-s*n)*U; ut = r**(1-s*n)*V
err = sp.diff(ur, r); ett = (ur+sp.diff(ut, th))/r
ert = sp.Rational(1,2)*(sp.diff(ur,th)/r + sp.diff(ut,r) - ut/r)
# shapes (factor r^{-sn})
erh = sp.simplify(err*r**(s*n)); eth = sp.simplify(ett*r**(s*n)); erh2 = sp.simplify(ert*r**(s*n))
# deviatoric strain shapes
d_rr = erh; d_tt = eth
e2 = sp.Rational(2,3)*(d_rr**2 + d_tt**2 + 2*erh2**2)
# constitutive (B=1): s_ij = (2/3) sigma_e^{1-n} eps_ij ; sigma_e^{1-n} = (sqrt(e2))^{(1-n)/n}
K = e2**((1-n)/(2*n))
s_rr = sp.Rational(2,3)*K*d_rr
s_tt = sp.Rational(2,3)*K*d_tt
s_rt = sp.Rational(2,3)*K*erh2
# stresses (shape): +pressure q = P(theta)
S_rr = s_rr + P; S_tt = s_tt + P; S_rt = s_rt
# equilibrium: r-part: r d(sig_rr)/dr + (sig_rr-sig_tt) + d(sig_rt)/dth =0 (sig=r^{-sn}S)
# => -sn S_rr + (S_rr - S_tt) + S_rt' = 0
eqr = sp.simplify(-s*n*S_rr + (S_rr - S_tt) + sp.diff(S_rt, th))
# theta-part: r d(sig_rt)/dr + 2 sig_rt + d(sig_tt)/dth =0 => -sn S_rt + 2S_rt + S_tt'=0
eqt = sp.simplify(-s*n*S_rt + 2*S_rt + sp.diff(S_tt, th))
# incompressibility already imposes her=... let's keep it as 3rd eq:
inc = sp.simplify((1-s*n)*U + U + sp.diff(V, th))  # err+ett shapes = 0
print("eqr=",eqr)
print("eqt=",eqt)
print("inc=",inc)