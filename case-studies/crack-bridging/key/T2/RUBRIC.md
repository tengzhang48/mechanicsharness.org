# Rubric — 100 points

Accept a beam boundary-value solution, an exact solution, or a demonstrably
equivalent numerical mechanics calculation. Code reuse and naming a source are
not required. Total moment on each arm is the requested convention.

For conceptual rows, give full points for a correct explicit argument or an
equivalent implemented calculation, half points for the right approach with a
missing necessary detail, and zero for an absent or contradictory treatment.
More specific subdivisions below take precedence.

| Mechanical description and adaptation | Points |
|---|---:|
| Uses elastic strip bending coupled to closing face forces, rather than an isolated traction-free crack with assigned toughness | 7 |
| Computes each arm's bending stiffness and adds opening compliances: $D_i=E_i h_i^3/12$, $C=1/D_1+1/D_2$, or the equivalent using total-width stiffnesses | 10 |
| Applies the actual softening curve locally to the relative gap and integrates its triangular area correctly | 7 |
| Enforces matching face displacement/slope at the front and retains nonzero intrinsic front work; e.g. $\delta''(0)=\sqrt{2C\Gamma_0}$ | 8 |
| Enforces zero rear shear for pure couples and obtains the moment from rear curvature with the correct width conversion | 8 |

These rows total **40 points**. For the stiffness row, award 5 points for correct
individual stiffnesses and 5 for their correct coupling. For the front row,
award 4 for displacement/slope matching and 4 for the intrinsic-work condition.
Equivalent energy or finite-element treatments need not write the displayed
equations.

| Numerical output | Target | Points |
|---|---:|---:|
| $M$, total on each arm | 0.212625787376721 N·m | 15 |
| $L$ | 10.8001023708973 mm | 15 |
| $\delta(L/2)$, full relative opening | 38.2188979260933 µm | 10 |

Numerical rows total **40 points**. For each, award its full weight for error
≤0.5%; 80% of its weight for error >0.5% to 2%; 50% for >2% to 5%;
20% for >5% to 15%; and zero above 15% or if missing. Use unrounded values
from [reference_results.json](reference_results.json). Fractional points are
allowed. Correctly converted SI outputs qualify.

Award **12 points for verification**: full credit for an executed, correctly
interpreted energy balance, independent exact/numerical route, constant-traction
or zero-front-energy limit, or meaningful convergence check. Give 6 for only
boundary/ODE residual checks, 3 for an unexecuted but relevant proposed check,
and 0 for none or for checks inconsistent with the answer.

Award **8 points for scope and event recognition**: 4 for using the prescribed
elastic, small-rotation regime up to the first complete connector failure, and
4 for recognizing that the zone/profile depends on both arms and geometry,
whereas the curve area alone determines neither length nor profile. Correct
formulation may demonstrate this implicitly; no extra essay is required.
Give 2 per part for an incomplete but consistent discussion, otherwise 0.

Correct numbers without method receive the numerical credit, not unshown
model or verification credit. Preserve conceptual credit for a sound adaptation
with arithmetic errors. In particular, identifying the energy area correctly
does not by itself earn the boundary-condition or zone-length modeling points.
