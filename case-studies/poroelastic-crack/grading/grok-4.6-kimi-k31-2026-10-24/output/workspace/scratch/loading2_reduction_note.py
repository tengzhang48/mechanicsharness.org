#!/usr/bin/env python3
"""Document the loading-2 long-time reduction error used by both submissions.

Both A and B replace the coupled drained problem
    sigma_n = 0, p = 1 on the faces; fields vanish at infinity
by a homogeneous drained elastic Griffith crack whose face load is the
Biot effective traction alpha.

Equilibrium is on TOTAL stress: sigma_ij,j = 0.  The Biot effective
stress bar_sigma_ij = sigma_ij + alpha p delta_ij then satisfies
    bar_sigma_ij,j = alpha p,j.
That body force vanishes only if p is spatially uniform.  Loading 2
holds p = 1 on a finite crack and p -> 0 at infinity, so p is not
uniform and the reduction is not the coupled problem.

The paper's long-time back-stress is the plane-strain modulus
    eta = alpha (1-2 nu) / (2 (1-nu)),
not alpha.  The two formulas below are identical; neither equals the
effective-traction substitution.
"""

import math

NU = 0.2
NU_U = 0.4
ALPHA = 0.89

eta = ALPHA * (1 - 2 * NU) / (2 * (1 - NU))

k2_inf_correct = -eta
v2_inf_eta = -eta * math.pi * (1 - NU)
v2_inf_alpha = -ALPHA * math.pi * (1 - 2 * NU) / 2.0

k2_inf_reduced = -ALPHA
v2_inf_reduced = -ALPHA * math.pi * (1 - NU)

print("eta =", eta)
print("correct K2_inf =", k2_inf_correct)
print("correct V2_inf (eta form)   =", v2_inf_eta)
print("correct V2_inf (alpha form) =", v2_inf_alpha)
print("eta-form matches alpha-form:", math.isclose(v2_inf_eta, v2_inf_alpha))
print("reduced  K2_inf =", k2_inf_reduced)
print("reduced  V2_inf =", v2_inf_reduced)
print("K ratio reduced/correct =", k2_inf_reduced / k2_inf_correct)
print("V ratio reduced/correct =", v2_inf_reduced / v2_inf_alpha)
print("combined K correct (load1+load2) =", 1.0 + k2_inf_correct)
print("combined K reduced  (load1+load2) =", 1.0 + k2_inf_reduced)
