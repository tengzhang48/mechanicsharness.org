import sympy as sp
s = sp.symbols('s')   # Fourier variable along x1 (careful: s also used as disp)
xi = sp.symbols('xi', real=True)
a, b = sp.symbols('a b', positive=True)

# ======================================================================
# Wiener-Hopf, semi-infinite Mode-I crack in strip period H on foundation.
# ======================================================================
# Use the anti-plane analogy first to fix the STRUCTURE, then map.
# Mode I in the strip: work with the strip [0,H], crack x2=0.
# Fourier transform in x1 (xi). For each xi, solve the x2-BVP on [0,H]
# to get the "face compliance" relation between opening transform
# Wtilde(xi) (on x1<0) and the traction transform on the crack plane.
#
# Dispersion: s^2 satisfies (derive3):
#   det = (C q^2 - C s^2 + kh)(C nu q^2 - C nu s^2 - C q^2 + C s^2 - 2kh)/2
# with q the x2-wavenumber. In the strip, x2-modes are q=pi n/H.
# The relevant branch for symmetric (mode I) opening is the one that
# connects continuously to the foundation pole at xi where s=0 -> q such
# that mode is marginal.
#
# Do the WH in the x1-Fourier variable xi directly for the SINGLE crack
# (H infinite), where modes in x2 are continuous, to get K for H->inf.
print("single crack WH")