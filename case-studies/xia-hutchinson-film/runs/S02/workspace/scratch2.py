import numpy as np
from scipy.integrate import quad

# ======================================================================
# CASE A: periodic semi-infinite cracks, period H, in film on foundation.
# Mode I. Use the antiplane/plane dislocation density method.
# ======================================================================
# Far behind the tips the strip is split; the film slides. Consider the
# opening profile w(x1) = u2(x1,0+) - u2(x1,0-) for x1<0 (one face pair).
# The crack faces carry no traction (x1<0). The film ahead is intact.
#
# Reduce to an integral equation on the crack line using the periodic
# climb-dislocation kernel. For a periodic array of climb dislocations
# of period H on x2=0, the normal stress sigma22 on the line is
#   sigma22(x1) = (Ebar/2) * PV ∫ B(x') K(x1-x') dx'
# where B = -d/dx'(opening) and the periodic kernel
#   K(s) = (pi/H) coth(pi s / H).
# (This is exact for the isolated crack plane pair interaction in an
#  infinite plane; it accounts for the neighbor cracks.)
# The foundation adds a body term; but far behind, the foundation is
# what stops the sliding. We treat the strip as: intact region x1>0 is
# elastic + foundation; the opening is determined by matching.
#
# THIS is getting unwieldy analytically. Do it NUMERICALLY: solve the
# strip problem by Fourier series in x2 (period H) and integration in x1.
print("switch to numerical BVP")