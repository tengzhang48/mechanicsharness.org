| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: antiplane compatibility, equilibrium, nonlinear shear law | 8 / 8 | Correct governing formulation: engineering shear is ∇w, equilibrium is div τ = 0, and the inverse power law gives the stated nonlinear equation. The implementation error is assessed under flux normalization below. |
| Model: admissible radial behavior | 8 / 8 | Correct exponent β = 1/(n+1) = 2/9, justified by finite, radius-independent nonzero J. |
| Model: odd displacement, traction-free faces, positive slip | 6 / 6 | The displacement is explicitly extended oddly, f′(±π) vanishes, and the upper-face slip is positive. The mistaken extension of its derivative causes the normalization error; it does not change these three stated conditions. |
| Model: magnitude fixed from J | 4 / 8 | Correct work density W = n/(n+1) τ·γ and correct flux expression, but one derivative-parity error corrupts the full-circle integral and hence A. This meets the rubric's one-normalization-error provision. |
| Numerical: Δw at 0.0200 mm | 4 / 8 | Reports +0.00183 mm versus +0.001318126778 mm; relative error 38.833%. Correct slip formula evaluated with the single erroneous amplitude earns half credit. |
| Numerical: Δw at 0.0800 mm | 4 / 8 | Reports +0.00249 mm versus +0.001793693738 mm; relative error 38.820%. Same upstream error; half credit. |
| Numerical: signed τ_rz at 0.0200 mm, 60° | 3 / 6 | Reports +42.4 MPa versus +38.58659686 MPa; relative error 9.883%. Correct component definition and sign, with the same amplitude error; half credit. |
| Numerical: signed τ_θz at 0.0200 mm, 60° | 3 / 6 | Reports +151.7 MPa versus +138.0106166 MPa; relative error 9.919%. Correct component definition and sign, with the same amplitude error; half credit. |
| Numerical: engineering-shear contour radius at 0° | 2 / 4 | Reports 0.07167 mm versus 0.04689442305 mm; relative error 52.833%. Correct contour formula with the same amplitude error; half credit. |
| Numerical: engineering-shear contour radius at 90° | 2 / 4 | Reports 0.03831 mm versus 0.02506612349 mm; relative error 52.836%. Same upstream error; half credit. |
| Numerical: engineering-shear contour radius at 180° | 2 / 4 | Reports 0.02048 mm versus 0.01339840659 mm; relative error 52.854%. Correct nonzero face limit and contour formula, with the same amplitude error; half credit. |
| Completion: seven identifiable numerical outputs | 5 / 5 | All seven requested quantities are present and identifiable. |
| Conventions: units, engineering shear, face limit | 5 / 5 | Correct mm/MPa units (2), engineering shear ∇w rather than tensor shear (2), and explicit upper-face limit (1). |
| Verification: genuine field/boundary requirement | 5 / 5 | The angular ODE's face derivative tests the traction-free condition. This check remains valid despite the flux-normalization error. |
| Verification: numerical result | 5 / 5 | The rerun reproduces f′(π) = 1.416×10⁻¹⁶ and near-face traction magnitude 4.153×10⁻⁸ MPa at r = 1 mm, consistent with the reported residual. |
| Verification: interpretation of agreement | 5 / 5 | Correctly interprets the residual as satisfying zero face traction. One sound boundary check earns all 15 points; the defective energy check supplies no additional credit. |
| Missing fracture physics | 5 / 5 | Correctly states that initiation requires a fracture/toughness criterion, such as a critical J, separation, or strain. |

**Total: 76 / 100.** Category totals: model 26/30; numerical correctness 20/40; completion and conventions 10/10; quantitative verification 15/15; missing fracture physics 5/5.

Scientific errors:

- In `submission/solve.py`, `uval(th)` returns `sgn(th)*sol.sol(abs(th))[1]`. Since f is odd, f′ must be even: the sign factor belongs in `fval`, but not in `uval`. At θ = −60°, a finite difference of `fval` gives f′ ≈ +0.7636947971, while `uval` supplies −0.7636947972. Thus the lower-half-plane shear used for the energy calculation is incompatible with the stated displacement.
- This single error cancels a necessary stress/displacement cross contribution in the full-circle flux. The reported coefficient J(A = 1) = 624.1363418 should be 953.8687488, and A = 0.001287364060 should be 0.000925609062 in the chosen normalization. The correctly evaluated flux for the submitted amplitude is 0.1833962264 N/mm, rather than the supplied 0.120 N/mm.
- Consequently, all requested slips, stresses, and contour radii are too large. Correcting only the derivative extension in a scratch copy restores all seven reference values to within 3.6×10⁻¹⁵ relative error. This establishes the single upstream error required for consequential-error credit; corrected results are diagnostic evidence, not credited submission outputs.

Unsupported claims:

- The statement that all seven reported numbers are within the 2% tolerance is false; both the printed answer and the unmodified code outputs exceed tolerance for every requested quantity.
- The reported energy checks do not establish that the compatible physical field carries the supplied J. Radius independence follows from the power scaling even with the erroneous derivative extension, and recovering the input J using the same erroneous coefficient does not validate that coefficient.

Presentation problems:

- The “Assumptions changed” section mostly repeats supplied assumptions rather than identifying changes.
- The “Missing / uncertain” discussion introduces unnecessary uncertainty about whether J fixes the local amplitude, despite the answer's correct governing argument that it does under the specified leading-field assumptions.

Tolerance: **0 of 7 requested numbers are within the stated 2% relative tolerance**; errors from the printed answer are 38.833%, 38.820%, 9.883%, 9.919%, 52.833%, 52.836%, and 52.854%, respectively.

Reference assessment: no reference error found. Its rerun reproduces the recorded values exactly, its two independent routes agree within 4.12×10⁻¹⁶ relative error, and its direct flux checks recover 0.120 N/mm. The one-line diagnostic repair independently reproduces its requested values, so no alternative corrected-reference grade is warranted. Rerun and diagnostic evidence is saved in `scratch/submission_solve_output.txt`, `scratch/submission_verify_output.txt`, `scratch/corrected_parity_output.txt`, `scratch/reference_rerun.json`, and `scratch/grading_checks.json`.
