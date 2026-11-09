# Grade

| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: compatibility, equilibrium, and constitutive law | 8 / 8 | Uses engineering shear γ = ∇w, div τ = 0, and the correct inverse nonlinear law, giving p = 9/7. |
| Model: admissible radial behavior | 8 / 8 | The correct a = 2/9 gives finite, nonzero, radius-independent J. The closed angular solution satisfies the nonlinear eigenproblem and is checked by independent ODE integration. |
| Model: symmetry, face traction, and positive slip | 6 / 6 | Odd F, even F′, F′(±π) = 0, and positive amplitude implement all three conventions (2 points each). |
| Model: magnitude fixed by J | 8 / 8 | Uses W = nτγ/(n+1) and both terms of the full flux in solve.py. J(A=1) = 137.9253071683 gives A = 0.00416524077788 in its angular normalization. |
| Numerical: Δw at 0.0200 mm | 8 / 8 | Reports +0.0013181 mm; relative error 0.002032%. Rerun reproduces the key. |
| Numerical: Δw at 0.0800 mm | 8 / 8 | Reports +0.0017937 mm; relative error 0.000349%. Rerun reproduces the key. |
| Numerical: signed τ_rz | 6 / 6 | Reports +38.59 MPa; relative error 0.008819%, with the correct sign. |
| Numerical: signed τ_θz | 6 / 6 | Reports +138.01 MPa; relative error 0.000447%, with the correct sign. |
| Numerical: contour radius at 0° | 4 / 4 | Reports 0.046894 mm; relative error 0.000902%. |
| Numerical: contour radius at 90° | 4 / 4 | Reports 0.025066 mm; relative error 0.000493%. |
| Numerical: contour radius at 180° | 4 / 4 | Reports 0.013398 mm; relative error 0.003035%. |
| Completion: all seven outputs | 5 / 5 | All requested values are numerical and clearly identifiable. |
| Conventions: units | 2 / 2 | Correct mm for slips/radii and MPa for stresses. |
| Conventions: engineering shear | 2 / 2 | Uses γ = ∇w and its magnitude without a tensor-shear factor-of-two error. |
| Conventions: face limit | 1 / 1 | Explicitly identifies the upper-face limit from within the material. |
| Verification: genuine requirement | 5 / 5 | Direct dimensional circle-flux integration checks the loading and path independence. Independent angular ODE integration provides another sound check. |
| Verification: numerical result | 5 / 5 | Rerun gives J = 0.120 N/mm at r = 0.005, 0.020, 0.100, and 0.300 mm (maximum absolute error 5.6×10⁻¹⁷ N/mm); ODE F′(π) = −2.835×10⁻¹⁷. |
| Verification: interpretation | 5 / 5 | Explains agreement with the imposed flux, proportional angular shape, and traction-free face condition. |
| Missing fracture physics | 5 / 5 | Correctly states initiation is undetermined and identifies a missing fracture criterion, such as critical J. |

**Total: 100 / 100** (model 30, numerical 40, completion/conventions 10, verification 15, missing physics 5).

The reference rerun passes its independent-route and limiting-solution assertions; no reference error was found. Rerun logs and the individual tolerance calculations are in scratch/submission.stdout, scratch/reference.stdout, and scratch/audit.json.

(a) Scientific errors:

- The statement that p-Laplacian homogeneity alone forces a = 2/9 overstates the argument. Finite nonzero flux requires 1 + p(a−1) = 0, or the slit angular eigenproblem must select the exponent. The correct field and eigenproblem verification are present, so radial-behavior credit remains full.
- The auxiliary submission/debug.py, line 26, omits amplitude A from w_,x. Evaluating that formula at the fitted amplitude gives J ≈ 16.4130064 N/mm instead of 0.120. This error is absent from solve.py, which supplies the correct requested results and credited flux verification.

(b) Unsupported claims:

- The quoted “relative ~10⁻⁹” equilibrium residual has no specified length normalization: div τ has units MPa/mm, while the compared stress has units MPa. The rerun's largest sampled absolute residual is 1.222×10⁻⁶ MPa/mm. This ambiguity does not affect the independently verified flux and angular-ODE checks.

(c) Presentation problems:

- Minor: an even stress function φ is mentioned without definition, and the “Assumptions changed” heading describes assumptions already supplied by the problem. Neither obscures the requested outputs or loses rubric points.

Tolerance: All seven requested numbers, including both signed stresses, are within the stated 2% relative tolerance; the largest relative error in the main displayed answers is 0.008819%.
