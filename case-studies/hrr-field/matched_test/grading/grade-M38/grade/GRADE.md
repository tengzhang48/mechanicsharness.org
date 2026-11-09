| rubric item | points awarded / max | reason |
|---|---:|---|
| Antiplane compatibility, equilibrium, and nonlinear law | 8 / 8 | Uses γ = ∇w, div τ = 0, and the correct inverse power law in both the governing PDE and code. The parenthetical law transcription noted below does not affect the implemented formulation. |
| Admissible local radial behavior | 8 / 8 | Solves the nonlinear angular eigenproblem with the face condition; rerun gives s = 0.222222222222202, agreeing with 1/(n+1) = 2/9 to 2.01×10⁻¹⁴. This is the admissible finite, nonzero-J exponent. |
| Odd symmetry, traction-free faces, and positive slip | 6 / 6 | Implements odd f, even f′, f′(π) = 0, and positive amplitude/upper-face displacement: 2 points for each requirement. |
| Magnitude fixed by J and correct stored work | 8 / 8 | Uses W = n/(n+1) τ·γ and both terms of the specified flux. Reproduces I_ang = 2.019595277308 and A = 9.256090617515×10⁻⁴ mm^(1−s) from J = 0.120 N/mm. |
| Δw at 0.0200 mm | 8 / 8 | Reported +0.0013181 mm versus +0.001318126778 mm; relative error 0.002032%, within 2%. |
| Δw at 0.0800 mm | 8 / 8 | Reported +0.0017937 mm versus +0.001793693738 mm; relative error 0.000349%, within 2%. |
| Signed τ_rz at 0.0200 mm, 60° | 6 / 6 | Reported +38.59 MPa versus +38.58659686 MPa; relative error 0.008819%, within 2%, with correct sign. |
| Signed τ_θz at 0.0200 mm, 60° | 6 / 6 | Reported +138.01 MPa versus +138.0106166 MPa; relative error 0.000447%, within 2%, with correct sign. |
| Engineering-shear contour radius at 0° | 4 / 4 | Reported 0.04689 mm versus 0.04689442305 mm; relative error 0.009432%, within 2%. |
| Engineering-shear contour radius at 90° | 4 / 4 | Reported 0.02507 mm versus 0.02506612349 mm; relative error 0.015465%, within 2%. |
| Engineering-shear contour radius at upper face, 180° | 4 / 4 | Reported 0.01340 mm versus 0.01339840659 mm; relative error 0.011893%, within 2%. Uses the nonzero limiting radial shear at the traction-free face. |
| Presence and identification of seven outputs | 5 / 5 | All seven requested numerical outputs are explicitly identified. |
| Units, engineering shear, and face-limit convention | 5 / 5 | Correct mm/MPa units (2), γ = ∇w as engineering shear (2), and explicit face-limit evaluation (1). |
| Quantitative verification | 15 / 15 | Genuine checks (5), numerical results (5), and agreement interpreted against the required conditions (5). Rerun circle integrals at r = 0.02 and 0.15 mm return 0.120000000 N/mm; the face traction ratio is 7.492×10⁻¹⁶; the n = 1 slip is 2.256758334×10⁻⁴ mm, matching the independent K-field formula. These suffice without crediting the unsupported finite-difference claim. |
| Missing fracture physics | 5 / 5 | Correctly states that initiation requires an additional criterion, identifying critical J, a finite-distance failure criterion, or a cohesive rule. |

**Total: 100 / 100.** The ancillary issues below do not invalidate the implemented model, requested outputs, or the sound quantitative checks that satisfy the rubric.

**(a) Scientific errors**

- ANSWER.md line 23 incorrectly restates the vector law as γ = γ₀(τ/τ₀)^n τ/τ₀. With the last factor denoting the stress vector divided by τ₀, the exponent must be n−1; alternatively, use exponent n with the unit stress direction τ/|τ|. The scalar inverse relation, PDE, work density, and code use the correct law.
- Line 36 reverses the outward face-traction signs: at θ = ±π the outward normals are ±e_θ, so t_z = ±τ_θz, rather than ∓τ_θz. This leaves the correctly imposed zero-traction condition and requested signed stresses unaffected.
- Lines 79–81 give an incorrect common load sensitivity. At fixed positions, A and slip scale as J^(7/9), stresses as J^(2/9), and the fixed-|γ| contour radii as J. They do not all change by the same relative amount as J.

**(b) Unsupported claims**

- The finite-difference divergence and curl residuals of order 3×10⁻¹⁰ (lines 69–70) are not calculated by the submitted code, and no evaluation points, step sizes, or residual output are supplied. They receive no verification credit.
- The assertion that real material data would deviate below γ ≈ 0.002 (lines 82–83) has no supporting data. γ₀ is a constitutive scaling parameter here, not an established elastic threshold.

**(c) Presentation problems**

- The linear-limit statement ends with the unfinished expression “K = 2J/...”. The complete relation is K = √(2μJ), which the preceding J relation and code correctly imply.
- The claimed ≤2.4×10⁻⁵ relative strain residual at the “reported radii” needs a rounding qualification. It holds for the six-decimal radii printed by solve.py (maximum 2.3603×10⁻⁵), while the coarser ANSWER.md table gives a maximum of 1.2027×10⁻⁴. Both are comfortably within the requested tolerance.

The submitted code and reference code both ran successfully in scratch/. The reference rerun reproduces the recorded results, and its two independent routes agree to a maximum relative discrepancy of 4.12×10⁻¹⁶. No reference error was found; a separate corrected-reference score is unnecessary.

**Tolerance: all seven requested numbers, including both stress signs, are within the stated 2% relative tolerance; the largest relative error is 0.015465%.**
