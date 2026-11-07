import sympy as sp

# ======================================================================
# Direct 1D reduction far behind (fully separated wake).
# Validate the two dispersion branches and identify ell.
# ======================================================================
# Plane stress, isotropic. Consider a state that varies only in x1:
# u1 = u(x1), u2 = 0  (uniaxial along x1)
# eps11 = u', others 0.
# sigma11 = Ebar[(1-nu)u' + nu u'] = Ebar u'
# sigma22 = Ebar[0 + nu u'] = Ebar nu u'
# equilibrium x1: h sigma11' - k u = 0 -> h Ebar u'' - k u = 0
#   -> decay length ell_1 = sqrt(h Ebar/k)
# This matches dispersion branch1 with q=0: s^2 = k/(h Ebar). GOOD.
Ebar = sp.symbols('Ebar', positive=True)
print("branch1 (u1 only):  ell = sqrt(h*Ebar/k)")

# Other 1D state: u2 = w(x1), u1=0
# eps22 = w', sigma22 = Ebar w', sigma12=0
# equilibrium x2: h sigma22' - k w = 0 -> h Ebar w'' - k w = 0
#   -> same ell = sqrt(h Ebar/k).  (mode I uniaxial)
print("branch2 (u2 only):  also ell = sqrt(h*Ebar/k)")

# So both diagonal components decay with the SAME length when the
# coupled term vanishes. But dispersion had two branches. The coupled
# (nu) branch is the one relevant to in-plane shearing both directions.
# Recompute branch2: s^2 = q^2 + 2k/(h Ebar (1-nu)). At q=0:
k,h,nu,E=sp.symbols('k h nu E',positive=True)
Ell_b2 = sp.sqrt(h*Ebar*(1-nu)/(2*k))
print("branch2 length:", Ell_b2)
# For nu=0.3, branch2 decays FASTER than branch1.
print("branch2/branch1 =", sp.sqrt((1-nu)/2))
print(" val nu=0.3:", float(sp.sqrt((1-0.3)/2)))