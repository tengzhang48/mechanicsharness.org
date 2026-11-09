| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: antiplane compatibility, equilibrium, and nonlinear shear law | 8 / 8 | The stress potential enforces equilibrium; the angular equation enforces compatibility of the supplied engineering-shear law. This is an equivalent valid formulation. |
| Model: admissible radial behavior | 8 / 8 | Correctly derives stress proportional to r^(-1/(n+1)), shear proportional to r^(-n/(n+1)), and displacement proportional to r^(1/(n+1)) from finite nonzero J. |
| Model: odd symmetry, traction-free faces, positive slip | 6 / 6 | Implements g′(0)=0, g(π)=0, and g(0)=-1. The resulting displacement is odd and positive on the upper face; the rerun gives face traction -3.201×10^-15 MPa. Award 2 for each convention. |
| Model: magnitude fixed by J | 8 / 8 | Code uses W=n/(n+1) τ·γ and the complete flux W cosθ−τ_rz w_,x. It computes C_J=2.84661855×10^-10 and a=82.53416679 in its stress-potential normalization from J=0.120 N/mm. |
| Numerical: Δw at 0.0200 mm | 8 / 8 | Reported +0.001318 mm; key +0.001318126778 mm. Relative error 0.009618%. |
| Numerical: Δw at 0.0800 mm | 8 / 8 | Reported +0.001794 mm; key +0.001793693738 mm. Relative error 0.017074%. |
| Numerical: signed τ_rz at 0.0200 mm, 60° | 6 / 6 | Reported +38.59 MPa; key +38.58659686 MPa. Relative error 0.008819%; sign correct. |
| Numerical: signed τ_θz at 0.0200 mm, 60° | 6 / 6 | Reported +138.01 MPa; key +138.0106166 MPa. Relative error 0.000447%; sign correct. |
| Numerical: engineering-shear contour radius at 0° | 4 / 4 | Reported 0.0469 mm; key 0.04689442305 mm. Relative error 0.011893%. |
| Numerical: engineering-shear contour radius at 90° | 4 / 4 | Reported 0.0251 mm; key 0.02506612349 mm. Relative error 0.135149%. |
| Numerical: engineering-shear contour radius at upper face | 4 / 4 | Reported 0.0134 mm; key 0.01339840659 mm. Relative error 0.011893%. |
| Completion: all seven numerical outputs | 5 / 5 | All seven values are present, identifiable, and reproduced by the submitted code. |
| Conventions: mm and MPa units | 2 / 2 | All requested slips/radii use mm and stresses use MPa. |
| Conventions: engineering versus tensor shear | 2 / 2 | Uses γ=∇w and γ₀=0.002 directly; no erroneous factor of two. |
| Conventions: face limit | 1 / 1 | Explicitly identifies the upper-face limit and evaluates the angular field at π from within the material. |
| Verification: genuine field or flux requirement | 5 / 5 | Recomputes the full J integral using Cartesian finite differences for w_,x on three circular arcs. Boundary and linear-limit checks supply additional valid verification. |
| Verification: numerical result or residual | 5 / 5 | Unchanged-code rerun gives J=0.120000 N/mm at r=0.020, 0.080, and 0.150 mm; g(π)=2.244×10^-17 and the linear amplitude check differs by 1.33×10^-15 relatively. |
| Verification: interpretation | 5 / 5 | Compares the calculated flux with the required 0.120 N/mm, the face residual with zero, and the linear result with J=K²/(2μ). These are meaningful agreement checks. |
| Missing fracture physics | 5 / 5 | Correctly states that initiation needs a separate fracture criterion, such as critical J or a critical separation, absent from the stationary field data. |

**Total: 100 / 100.** Category subtotals: model 30/30; numerical correctness 40/40; completion/conventions 10/10; verification 15/15; missing fracture physics 5/5.

The submitted solver was copied unchanged to `scratch/submission_solve.py` and completed successfully; its output is in `scratch/submission_run.txt`. The reference rerun passed all assertions and reproduced the recorded key, including its independent phase solution. No reference error was found, so an alternate corrected-reference grade is unnecessary. Individual comparisons use the values printed in ANSWER.md, including their rounding; details are in `scratch/number_comparison.json`.

**(a) Scientific errors**

- No error in the governing stationary field, work density, amplitude, or requested outputs was found.

**(b) Unsupported claims**

- The claim of independently determining the radial eigenvalue is overstated: the shooting code varies m while keeping α fixed at the already derived value 2/9. A general stress-potential compatibility equation requires α=1+n(m−1). This checks consistency with the imposed scaling rather than independently deriving that scaling. The correct finite-J argument and the reproduced energy/boundary/linear checks still earn full model and verification credit.

**(c) Presentation problems**

- The auxiliary finite-difference residuals labeled “relative” retain units of inverse mm: the code divides stress/strain gradients by stress/strain magnitudes without a length factor. Multiplying by r=0.03 mm gives dimensionless residuals approximately 2.27×10^-8 and 4.60×10^-9. This labeling issue does not affect the requested units or the full-credit quantitative checks.

Tolerance: **All seven requested numbers are within the stated 2% relative tolerance, with correct signs; the largest reported-value error is 0.135149%.**
