import sympy as sp
import time

# Derive the plane-strain incompressible power-law compatibility equation for
# the Airy stress function Phi = r^{2-s} f(theta).
# eps_ij = C g^{n-1} (dev) with g = sigma_e * r^s.
s, n = sp.symbols('s n')
th = sp.symbols('theta')
f = sp.Function('f')

Sig_rr = sp.diff(f(th), th, 2) + (2 - s)*f(th)
Sig_tt = (2 - s)*(1 - s)*f(th)
Sig_rt = (s - 1)*sp.diff(f(th), th)

# 3rd dimension: plane strain incompressible -> szz = 0 (deviatoric), sig_zz=1/2(srr+stt)
# deviatoric shapes:
D = sp.simplify(Sig_rr - Sig_tt)
# mean: (Sig_rr+Sig_tt+Sig_zz)/3 ; Sig_zz halves -> use deviatoric:
srr = sp.Rational(1,2)*D
stt = -sp.Rational(1,2)*D
srt = Sig_rt
ge2 = sp.Rational(3,2)*(srr**2 + stt**2 + 2*srt**2)
ge = sp.sqrt(sp.simplify(ge2))
# strains
C = sp.symbols('C', positive=True)
err = C*ge**(n-1)*srr
ett = C*ge**(n-1)*stt
ert = C*ge**(n-1)*srt

# 2D compatibility (strain from a displacement field):
#  r^2 d2(ett)/dr^2 + 2 r d(ett)/dr - r d(err)/dr + d2(err)/dth^2 - 2*d(ert)/dth - ett ... 
# Standard polar compatibility:
#  (1/r^2) d2( err... ) -> use the well-known form:
#  d2(ett)/dr^2*(r^2) + 2r d(ett)/dr - r d(err)/dr + d2(err)/dth^2
#      - 2 d(ert)/dth - r d2( (r ert) )... 
# Use the displacement-based compatibility instead (safer): define u_r,u_theta.

Ur = sp.Function('Ur')(th); Ut = sp.Function('Ut')(th)
r = sp.symbols('r', positive=True)
ur = r**(1 - s*n)*Ur
ut = r**(1 - s*n)*Ut
err_k = sp.diff(ur, r)
ett_k = (ur + sp.diff(ut, th))/r
ert_k = sp.Rational(1,2)*(sp.diff(ur, th)/r + sp.diff(ut, r) - ut/r)
# strain shapes: err_k = r^{-sn} * (1-sn)Ur ; ett_k = r^{-sn}(Ur+Ut') ; ert_k = r^{-sn}/2 (Ur'-sn Ut)
err_h = sp.simplify(err_k * r**(s*n))
ett_h = sp.simplify(ett_k * r**(s*n))
ert_h = sp.simplify(ert_k * r**(s*n))
print("err_h", err_h)
print("ett_h", ett_h)
print("ert_h", ert_h)
# constitutive inverts: given the strain field, need s_ij = (2/3) sigma_e^{1-n} eps_ij /C'
# with sigma_e = (eps_e / C')^{1/n}, eps_e^2=(2/3)(err^2+ett^2+2ert^2). So
# s_ij = (2/3)(eps_e^{2})^{(1-n)/(2n)} eps_ij / C'^{...}; we only need ratios.
epse2 = sp.Rational(2,3)*(err_h**2 + ett_h**2 + 2*ert_h**2)
H = epse2**sp.Rational(1-n, 2*n)   # = sigma_e^{1-n} * const
srr_h = sp.Rational(2,3)*H*err_h
stt_h = sp.Rational(2,3)*H*ett_h
srt_h = sp.Rational(2,3)*H*ert_h
sig_rr = srr_h + sp.Symbol('P')     # + pressure (hydrostatic) term
sig_tt = stt_h + sp.Symbol('P')
# equilibrium: two equations; eliminate P by cross-differentiating
r2 = sp.symbols('r')
# equilibrium in terms of shapes: sigma_ij = r^{-sn} * (shape) + P(r,th) delta
# d(sig_rr)/dr + (sig_rr-sig_tt)/r + d(sig_rt)/dth /r = 0
# We'll build with full r dependence:
sig_rr_full = r**(-s*n)*srr_h + sp.Function('P')(r, th)
print("built; now this route needs P elimination -- stopping (use displacement eqs instead)")