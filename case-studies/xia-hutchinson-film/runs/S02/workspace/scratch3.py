import numpy as np
from scipy.linalg import solve
from scipy.optimize import brentq

# ======================================================================
# CASE A numerical: periodic semi-infinite cracks, strip period H.
# Solve the film-foundation problem on a finite domain x1 in [-Xf, Xa],
# x2 in [-H/2, H/2] (one period), cracks at x2=±H/2 for x1 in [-Xf,0].
# Unknowns: (u1,u2) nodal. Boundary: periodic in x2.
# Use symmetric (mode I) reduction? Mixed mode may appear because cracks
# at x2=±H/2 and equibiaxial loading; but a periodic crack array under
# equibiaxial far field gives pure mode I at the tips by symmetry
# (reflection x1-> -x1? not a symmetry of semi-infinite). Let's just
# solve full and read K from near-tip field.
# ======================================================================
# Strategy: use superposition. Far field u_ahead=0 (sigma=s0 I).
# Crack faces x1<0, x2=±H/2 are traction-free.
# Solve: h sig_ab,b - k u_a = 0 in film, with sigma=s0 I + Dsig.
# Traction-free on crack faces. Periodic in x2, bounded ahead/behind.
print("Use Fourier-in-x2 / Wiener-Hopf instead.")