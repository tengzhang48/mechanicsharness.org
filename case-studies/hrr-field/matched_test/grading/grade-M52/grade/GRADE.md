| rubric item | points awarded / max | reason |
|---|---:|---|
| Antiplane compatibility, equilibrium, and nonlinear shear law | 8 / 8 | Uses γ = ∇w, div τ = 0, and the supplied pure power law. The stress function enforces equilibrium, and its angular ODE correctly enforces compatibility. |
| Admissible local radial behavior | 8 / 8 | Derives stress exponent s = 1/(n+1) = 2/9 from radius-independent, finite nonzero J. The resulting displacement exponent is also 2/9; the angular boundary residual confirms admissibility. |
| Symmetry, face traction, and positive slip | 6 / 6 | Φ is even and Φ′ is odd, giving odd w; Φ(π) = 0 enforces zero face traction; Φ(0) = −1 gives positive upper-face displacement. All three conditions are implemented correctly (2 points each). |
| Magnitude from J and correct stored work | 8 / 8 | Uses W = n/(n+1) τ·γ and both terms of the prescribed flux. The stress normalization gives I = 0.65179982 and B = 82.53416679 MPa·mm^(2/9), reproduced by the code. |
| Δw at 0.0200 mm | 8 / 8 | Main reported value +0.001318 mm versus +0.001318126778 mm: relative error 0.009618%. Parenthesized value and code also pass. |
| Δw at 0.0800 mm | 8 / 8 | Main reported value +0.001794 mm versus +0.001793693738 mm: relative error 0.017074%. Parenthesized value and code also pass. |
| Signed τ_rz at 0.0200 mm, 60° | 6 / 6 | +38.59 MPa versus +38.58659686 MPa: relative error 0.008819%, with correct sign. |
| Signed τ_θz at 0.0200 mm, 60° | 6 / 6 | +138.0 MPa versus +138.0106166 MPa: relative error 0.007693%, with correct sign. |
| Engineering-shear contour radius at 0° | 4 / 4 | Main reported value 0.0469 mm versus 0.04689442305 mm: relative error 0.011893%. Parenthesized value and code also pass. |
| Engineering-shear contour radius at 90° | 4 / 4 | Main reported value 0.0251 mm versus 0.02506612349 mm: relative error 0.135149%. Parenthesized value and code also pass. |
| Engineering-shear contour radius at 180° | 4 / 4 | Main reported face-limit value 0.0134 mm versus 0.01339840659 mm: relative error 0.011893%. Parenthesized value and code also pass. |
| All seven outputs present and identifiable | 5 / 5 | All requested values are explicitly tabulated. |
| Units, engineering shear, and face limit | 5 / 5 | Correct mm/MPa units (2), engineering shear γ = ∇w with no factor-of-two error (2), and explicit upper-face limit (1). |
| Verification addresses a genuine requirement | 5 / 5 | Checks traction-free faces, dimensional J on several circles, and the independently known n = 1 solution. |
| Verification supplies numerical results | 5 / 5 | Rerun gives Φ(π) = 1.478×10⁻¹⁶, face traction −2.263×10⁻¹⁴ MPa, and J = 0.120000000 N/mm at r = 0.020, 0.080, and 0.200 mm. Linear-limit agreement is at floating-point precision. |
| Verification interprets agreement | 5 / 5 | Correctly identifies vanishing face traction, radius independence and agreement with the supplied J, and agreement with the exact linear field. |
| Recognition of missing fracture physics | 5 / 5 | States that initiation requires an additional growth criterion, such as critical J or a cohesive rule; the stationary field alone does not determine initiation. |

**Total: 100 / 100** (model 30/30; numbers 40/40; completion/conventions 10/10; verification 15/15; missing physics 5/5).

Both the submission and reference were rerun from copies in `scratch/`. The reference's independent routes and assertions passed; no reference error was found, so a separate corrected-reference grade is unnecessary. Rerun evidence is in [submission_run.txt](scratch/submission_run.txt), [reference_run.json](scratch/reference_run.json), and [grade_checks.json](scratch/grade_checks.json). The minor observations below do not meet any rubric condition for losing points.

(a) Scientific errors:

- The ancillary decimal ratio is incorrect: 4^(2/9) = 1.360790000174377, not 1.3617 “exactly.” The formula, radial scaling, and both requested slips are correct.

(b) Unsupported claims:

- The suggestion of possible few-percent outer-geometry corrections at r = 0.08 mm has no supplied evidence and conflicts with the problem's instruction that corrections are negligible at all requested positions. The actual calculation nevertheless uses the required leading field.
- The description of s* as a numerically reproduced shooting eigenvalue overstates the code's procedure. With the small residual obtained here, `solve()` assigns `sstar = s0`; it checks the face condition at the exponent derived from J rather than independently finding the exponent. This is still a valid admissibility check and radial-exponent argument.

(c) Presentation problems:

- None material: the seven outputs, signs, units, and face-limit convention are clear. The optional uncertainty discussion could be shorter and adhere more closely to the stipulated leading-field assumption.

Tolerance: **All seven requested numbers are within the stated 2% relative tolerance**, including the signed stresses; the largest error among the main displayed values is 0.135149%.
