| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: antiplane compatibility, equilibrium, nonlinear shear law | 8 / 8 | Correctly gives γ = ∇w, ∇·τ = 0, and the inverse supplied power law. The submitted code implements the correct angular equation despite the transcription error described below. |
| Model: admissible radial behavior | 8 / 8 | Derives the stress exponent and displacement exponent q = 1/(n+1) = 2/9 from finite, nonzero J; the shooting calculation confirms this exponent. |
| Model: symmetry, traction-free faces, positive slip | 6 / 6 | Implements odd F, even F′, zero face traction, and positive upper-face displacement: 2 points each. |
| Model: magnitude fixed by J | 8 / 8 | Uses W = n/(n+1) τ·γ and both terms in the prescribed flux. Rerun gives ℐ = 2.019595277308 and C = 9.256090617515×10⁻⁴, consistent with the reference. |
| Numerical: Δw at 0.0200 mm | 8 / 8 | +0.0013181 mm; relative error 0.00203150%, within 2%. |
| Numerical: Δw at 0.0800 mm | 8 / 8 | +0.0017937 mm; relative error 0.00034910%, within 2%. |
| Numerical: signed τ_rz at 0.0200 mm, 60° | 6 / 6 | +38.5866 MPa; relative error 0.00000813%, within 2%, with correct sign. |
| Numerical: signed τ_θz at 0.0200 mm, 60° | 6 / 6 | +138.0106 MPa; relative error 0.00001205%, within 2%, with correct sign. |
| Numerical: engineering-shear contour radius at 0° | 4 / 4 | 0.046894 mm; relative error 0.00090214%, within 2%. |
| Numerical: engineering-shear contour radius at 90° | 4 / 4 | 0.025066 mm; relative error 0.00049265%, within 2%. |
| Numerical: engineering-shear contour radius at 180° | 4 / 4 | 0.013398 mm; relative error 0.00303459%, within 2%. |
| Completion: seven identifiable numerical outputs | 5 / 5 | All seven are clearly tabulated. |
| Conventions: mm/MPa units | 2 / 2 | All requested outputs use the specified units. |
| Conventions: engineering versus tensor shear | 2 / 2 | Uses γ = ∇w and its magnitude, without a spurious factor of two. |
| Conventions: face limit | 1 / 1 | Evaluates the angular field at π from the material interval; the nonzero face shear magnitude is retained. |
| Verification: genuine requirement | 5 / 5 | Computes the dimensional energy flux directly from stress and displacement gradients on several circular contours. |
| Verification: numerical result | 5 / 5 | Rerun of verify.py gives J = 0.12000000 N/mm at R = 0.01, 0.02, 0.05, and 0.20 mm for n = 3.5. It also reproduces the stated equilibrium and constitutive residuals. |
| Verification: interpretation | 5 / 5 | Explicitly identifies agreement with the supplied J and independence of contour radius. This sound check earns all 15 verification points. |
| Missing fracture physics | 5 / 5 | States initiation is undetermined without a fracture criterion, such as J_c, and a relation between applied load and J. |

**Total: 100 / 100** (model 30, numerical 40, completion/conventions 10, verification 15, missing physics 5).

The reference solver was rerun successfully; its independent angular and phase routes agree. No reference error was found. Full model credit follows the rubric's allowance for any correct formulation: the correct governing PDE, stress formulas, and implemented solver establish the model despite the written angular-equation typo. Evidence is retained in `scratch/`.

**(a) Scientific errors**

- In ANSWER.md, the angular multiplier is written H = (q²F² + F′²)^(1/n − 1). It must be H = (q²F² + F′²)^((1/n − 1)/2), equivalently g^(1/n − 1). The missing factor of one-half is inconsistent with the stated stress law. All submitted solvers implement the correct expression, so it does not affect the requested outputs.

**(b) Unsupported claims**

- The claim that face stress τ_θ(π) is zero to 10⁻¹⁶ does not specify radius, units, or normalization. The reproduced ≈−6×10⁻¹⁷ residual is F′(π), not dimensional stress; solve.py's field gives τ_θ ≈−1.82×10⁻¹⁴ MPa at r = 0.0200 mm. Face traction is numerically negligible, but the stated stress precision is not established as written.
- The scripts do not implement the claimed independent comparison with J = K²/(2μ), or explicitly compare the n = 1 angular solution with 2 sin(θ/2). The quoted linear slip is nevertheless correct: independent evaluation gives 0.000552790639154 mm. The demonstrated contour-flux check already satisfies the verification rubric.

**(c) Presentation problems**

- The unexplained “ν →” in the constitutive-law sentence is a typographical distraction.

**Tolerance:** All seven requested signed numerical values are within the stated 2% relative tolerance; the largest relative error is 0.00303459% at the upper-face contour radius.
