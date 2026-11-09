| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: antiplane compatibility, equilibrium, and nonlinear law | 8 / 8 | Uses engineering shear γ = ∇w, ∇·τ = 0, and the correct inverse of the supplied pure power law; the angular equilibrium equation is correct. |
| Model: admissible radial behavior | 8 / 8 | Derives a finite, nonzero flux from J ∝ r^(1−p(n+1)), giving stress and displacement exponents p = s = 1/(n+1) = 2/9. Here p is the submission's stress exponent. |
| Model: symmetry, free faces, and positive slip | 6 / 6 | Implements odd displacement and positive f′(0); final.py integrates both halves separately and gives f(+π) = −f(−π) = 1.698433234023 and f′(±π) = −3.40×10⁻¹⁵. Upper-face displacement is positive. |
| Model: magnitude fixed by J and correct stored work | 8 / 8 | final.py and verify2.py use W = n/(n+1) τ·γ and the full flux term W cosθ − τ_rz w_,x. Rerunning gives I_J = 953.868749 and A = 9.256090617515×10⁻⁴, with the material coefficient included in I_J. |
| Numerical: Δw at 0.0200 mm | 8 / 8 | Reports +0.00131813 mm versus +0.00131812677776 mm; relative error 0.000244456%. |
| Numerical: Δw at 0.0800 mm | 8 / 8 | Reports +0.00179369 mm versus +0.00179369373814 mm; relative error 0.000208405%. |
| Numerical: τ_rz at 0.0200 mm, 60° | 6 / 6 | Reports +38.587 MPa versus +38.5865968623 MPa; relative error 0.001044761%. |
| Numerical: τ_θz at 0.0200 mm, 60° | 6 / 6 | Reports +138.011 MPa versus +138.010616632 MPa; relative error 0.000277781%. |
| Numerical: contour radius at 0° | 4 / 4 | Reports 0.0468944 mm versus 0.0468944230549 mm; relative error 0.000049164%. |
| Numerical: contour radius at 90° | 4 / 4 | Reports 0.0250661 mm versus 0.0250661234889 mm; relative error 0.000093708%. |
| Numerical: contour radius at upper face, 180° | 4 / 4 | Reports 0.0133984 mm versus 0.0133984065871 mm; relative error 0.000049164%. |
| Completion: all seven numerical outputs | 5 / 5 | All seven quantities are present and identifiable. Both requested stress signs are explicit. |
| Conventions: units, engineering shear, and face limit | 5 / 5 | Correct mm/MPa units (2), γ = ∇w rather than tensor shear (2), and explicit upper-face limit from the material (1). |
| Verification: genuine field or energy requirement | 5 / 5 | Checks flux at multiple radii, Cartesian equilibrium, the constitutive law, face traction, and odd symmetry. These are actual calculations in final.py and verify2.py. |
| Verification: quantitative results | 5 / 5 | Rerun verify2.py gives J = 0.1200000000 N/mm at r = 0.02, 0.08, 0.20, and 0.30 mm; equilibrium residual/max terms = 2.571×10⁻¹⁰, constitutive residual = 2.550×10⁻¹⁶, and face traction = −6.266×10⁻¹³ MPa. |
| Verification: interpretation | 5 / 5 | Correctly interprets the flux as radius-independent and equal to the supplied J, and the small residuals as agreement with equilibrium, the shear law, and the traction-free condition. |
| Missing fracture physics | 5 / 5 | States that initiation requires a fracture criterion, identifying critical J/toughness or a critical strain/process-zone criterion as missing. |

**Total: 100 / 100.** Full credit follows from the final answer and the corrected final.py/verify2.py implementations. No consequential-error allocation is needed.

**(a) Scientific errors**

- No scientific error was found in the final reported field or requested outputs. The reference rerun reproduces its recorded values exactly and passes its independent phase-solution and linear-limit assertions; no corrected-reference alternative grade is warranted.
- The earlier solve2.py evaluates its angular solution at negative angles outside the integrated interval instead of making an odd extension. Checking its initialization gives f(−π) = 447.04412555 rather than −1.69843323402, I_J = 9617.95452019 rather than 953.868749, and an incorrect amplitude. final.py and verify2.py correct this, so the draft error does not reduce credit for the final solution.
- The earlier verify.py uses the angular magnitude g in the inverse stress law where the physical magnitude A r^(s−1) g is required. Its stress and flux checks therefore use an incorrectly scaled stress field. verify2.py corrects that error.

**(b) Unsupported claims**

- The claimed numerical n = 1 solver run has no corresponding implementation or saved run in the submitted files. The stated analytic linear limit is correct; the reproduced flux and residual checks independently earn full verification credit.
- The asserted overall quadrature/ODE uncertainty of 10⁻¹⁰ is not established by a refinement or convergence study. Small residuals and solver tolerances alone do not supply that error bound.

**(c) Presentation problems**

- Earlier scripts remain alongside the corrected implementations, and solve2.py is cited as a first pass despite producing a conflicting normalization. Reproduction should use the identified final.py and verify2.py. This does not obscure the seven final outputs or warrant a rubric deduction.

**Tolerance:** All seven requested signed numbers are within the stated 2% relative tolerance; the largest relative error in the displayed answer is 0.001044761%.
