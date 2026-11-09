# Grade

| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: compatibility, equilibrium, nonlinear shear law | 8 / 8 | Uses γ = ∇w, div τ = 0, and the correct inverse law τ = C\|γ\|^(1/n−1)γ in the answer and code. The incorrect label for the p-Laplacian index does not alter these equations. |
| Model: admissible radial behavior | 8 / 8 | Correct finite, nonzero J scaling gives w ∝ r^(2/9), γ ∝ r^(−7/9), and τ ∝ r^(−2/9). The angular equation and exponent-perturbation check are consistent. |
| Model: odd symmetry, traction-free faces, positive slip | 6 / 6 | Odd f and even f′ implement symmetry (2); f′(π) is numerically zero (2); positive A produces the stipulated upper-face motion and positive slip (2). |
| Model: magnitude fixed by J | 8 / 8 | Correct W = (7/9)τ·γ and both terms of the flux integral give I = 2.0195952773 and A = 9.2560906175 × 10⁻⁴ mm^(7/9). |
| Numerical: Δw at 0.0200 mm | 8 / 8 | Reported +0.0013181 mm; relative error 0.00203%, within 2%. |
| Numerical: Δw at 0.0800 mm | 8 / 8 | Reported +0.0017937 mm; relative error 0.000349%, within 2%. |
| Numerical: signed τ_rz at 0.0200 mm, 60° | 6 / 6 | Reported +38.59 MPa; relative error 0.00882%, within 2%; sign correct. |
| Numerical: signed τ_θz at 0.0200 mm, 60° | 6 / 6 | Reported +138.0 MPa; relative error 0.00769%, within 2%; sign correct. |
| Numerical: contour radius at 0° | 4 / 4 | Reported 0.04689 mm; relative error 0.00943%, within 2%. |
| Numerical: contour radius at 90° | 4 / 4 | Reported 0.02507 mm; relative error 0.01547%, within 2%. |
| Numerical: contour radius at upper face, 180° | 4 / 4 | Reported 0.01340 mm; relative error 0.01189%, within 2%. |
| Completion: all seven outputs | 5 / 5 | All requested numerical quantities are present and identifiable. |
| Conventions: units, engineering shear, face limit | 5 / 5 | Correct mm/MPa units (2), γ = ∇w with engineering-shear magnitude (2), and explicit upper-face limit evaluated at θ = π (1). |
| Verification: genuine requirement | 5 / 5 | Checks traction-free faces, rectangular-contour energy flux, Cartesian equilibrium, and the known n = 1 limit. |
| Verification: quantitative result | 5 / 5 | Rerun gives normalized face residual 6.446 × 10⁻¹⁶ and J = 0.120000 N/mm on both rectangles. The n = 1 angular error is 3.71 × 10⁻¹³. |
| Verification: interpretation | 5 / 5 | Explicitly compares residuals with zero, flux with the imposed J, and the linear field with its known solution. These sound checks suffice despite the minor derivative-bound overstatement below. |
| Missing fracture physics | 5 / 5 | Correctly states that initiation requires an additional fracture criterion, giving critical J and cohesive/process-zone examples. |

**Total: 100 / 100** — model 30/30; numerical correctness 40/40; completion and conventions 10/10; verification 15/15; missing fracture physics 5/5.

The reference rerun passed its assertions: the angular-ODE and independent phase solutions agree, with maximum relative difference 4.12 × 10⁻¹⁶. No reference error was found; no alternative grade is needed. The submitted solver reproduces all seven reference values. Rerun outputs and the per-number comparison are retained in `scratch/`.

Scientific errors:

- The answer calls its governing equation a p-Laplacian with p = 1/n. Under the standard definition div(\|∇w\|^(p−2)∇w) = 0, the correct index is p = 1 + 1/n = 9/7. Its displayed shear law and implemented angular equation are nevertheless correct, so the rubric's governing-model requirements are met.
- “Faces open: Δw > 0” confuses tangential antiplane slip with normal opening. With only u_z nonzero, positive slip does not demonstrate normal separation. This incidental statement does not change the correctly implemented loading convention.

Unsupported claims:

- The stated finite-difference versus analytic γ_x agreement of 1 × 10⁻¹¹ is slightly overstated: the two rerun relative errors are 7.94 × 10⁻¹² and 1.73 × 10⁻¹¹. The boundary, energy, and limiting-solution checks remain sound and earn full verification credit.
- Calling small strain “marginal” at about 2% shear has no quantitative validity threshold supplied in the problem or established by the answer. The reported strain magnitude and its divergence as r → 0 are supported.

Presentation problems:

- The symbol p is reused for the nonlinear operator index and the stress singularity exponent, obscuring otherwise correct radial-scaling reasoning.
- The supplementary slip ratio is printed as 1.3606; 4^(2/9) = 1.36079000017. This is not one of the seven requested numbers and does not affect their credit.

Tolerance verdict: **All seven requested numbers are within the stated 2% relative tolerance, with correct signs; the largest reported-value error is 0.01547%.**
