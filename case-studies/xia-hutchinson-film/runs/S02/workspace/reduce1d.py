import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import sympy as sp

# ======================================================================
# Exact 1D boundary-layer model, CASE A.  (reduced model)
# ======================================================================
# Far ahead (x1->+inf): intact, sigma = s0 I, u=0.
# Behind (x1<0): crack faces traction free.
# Foundation ties the released displacement back to zero with length ell.
#
# Reduce with the strip's fundamental "relaxation" compliance. Model the
# film as a 1D bar along x1 of effective modulus (plane stress uniaxial):
#   e11 = u1' ; sig11 = s0 + C e11   (C=E/(1-nu^2)? for uniaxial it is E)
# For the film between two cracks (period H), the average x1-strain
# relations combine with the opening.
#
# The cleanest self-consistent reduced model (standard wide-slab):
# The change of potential energy when material point relaxes from intact
# to fully split, while restrained by foundation over length ell, is
#   per period H per unit advance: U = (1/2) s0^2 / Ebar * H * L(H)
# with a dimensionless H-dependent factor from foundation. Then G = U/h.
#
# We now DERIVE this factor from the ODE by energy balance on the wake:
# behind the tip, the released displacement u1 (sliding) satisfies
#   C h u1'' - k u1 = 0   (C=Ebar),  u1(0)=0 (built-in at tip, intact
#   continuity), u1(-inf) -> free value where sig11=0 => Ebar u1' = -s0.
# Solve: u1(x1)= A e^{x1/ell}, x1<0, but requirement sig11=s0+Ebar u1'=0
# at x1->-inf is incompatible with exponential decay (u1'->0). So the
# sliding saturates only through the neighbor/interaction. 
# => the wake is NOT a single decaying mode; the release is provided by
#    the crack opening (mode I), not sliding. Let us instead use mode-I:
#    opening w(x1) with sig22 ahead = s0, w(0)=0, and behind the effective
#    stiffness of the film+foundation gives w'' = w/ell^2.
#    w(x1) = W e^{x1/ell}? again w'->0 => sig22=s0 not 0. Contradiction.
#
# RESOLUTION: The far-behind state is NOT stress-free in the sigma22
# sense; it is the state where the CRACK FACES carry zero traction, which
# for the film means the film relaxes the *x2* normal stress to zero at
# the faces but it is still connected ahead. The correct far-behind
# relaxed state is uniaxial tension in x1 with sig11=s0? No...
print("see next: use energy method with strip modes")