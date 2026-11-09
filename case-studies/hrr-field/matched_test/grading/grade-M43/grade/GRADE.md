# Grade

| rubric item | points awarded / max | reason |
|---|---:|---|
| Antiplane compatibility, equilibrium, and nonlinear shear law | 8 / 8 | States γ = ∇w, div τ = 0, and the correct inverse power law; implements a valid upper-half stress-function solution. The extension error is assessed under symmetry below. |
| Admissible local radial behavior | 8 / 8 | Derives the finite, nonzero-flux condition 1 − p(n+1) = 0 for stress exponent p, obtaining w ∝ r^(2/9) and τ ∝ r^(−2/9). |
| Odd displacement symmetry | 1 / 2 | Correctly states and implements odd w in `w_of`, but `tau_cart` and the circular-flux routines apply the wrong reflection signs to stresses and shear. Award partial credit for correct intent with an implementation transcription error. |
| Traction-free faces | 2 / 2 | Angular solution has τ_θz → 0 at θ → ±π; rerun face residual is 3.282 × 10⁻⁶ MPa at r = 0.05 mm. |
| Positive-slip loading convention | 2 / 2 | Positive amplitude gives positive upper-face displacement and positive relative face slip. |
| Magnitude fixed by J with correct work density and flux | 8 / 8 | Uses W = n τ·γ/(n+1) and both flux terms; computes J(A=1) = 2.84661855 × 10⁻¹⁰ and A = 82.53416679. The reflection error reverses both factors in τ_rz w_,x and leaves this integral and amplitude unchanged. |
| Δw at 0.0200 mm | 8 / 8 | Reports +0.0013181 mm versus +0.001318126778 mm; relative error 0.002032%. |
| Δw at 0.0800 mm | 8 / 8 | Reports +0.0017937 mm versus +0.001793693738 mm; relative error 0.000349%. |
| Signed τ_rz at 0.0200 mm, 60° | 6 / 6 | Reports +38.59 MPa versus +38.58659686 MPa; relative error 0.008819%, with correct sign. |
| Signed τ_θz at 0.0200 mm, 60° | 6 / 6 | Reports +138.01 MPa versus +138.0106166 MPa; relative error 0.000447%, with correct sign. |
| Engineering-shear contour radius at 0° | 4 / 4 | Reports 0.04689 mm versus 0.04689442305 mm; relative error 0.009432%. |
| Engineering-shear contour radius at 90° | 4 / 4 | Reports 0.02507 mm versus 0.02506612349 mm; relative error 0.015465%. |
| Engineering-shear contour radius at 180° | 4 / 4 | Reports 0.01340 mm versus 0.01339840659 mm; relative error 0.011893%. |
| All seven numerical outputs present | 5 / 5 | All requested quantities are clearly identified in the answer table. |
| Correct mm/MPa units | 2 / 2 | Slip and radii use mm; stresses use MPa. |
| Engineering versus tensor shear | 2 / 2 | Uses γ = ∇w and the engineering-shear magnitude 0.0100, without an erroneous factor of two. |
| Explicit face-limit convention | 1 / 1 | Labels the 180° result as the face limit and evaluates the nonzero limiting T(π). |
| Verification addresses a genuine requirement | 5 / 5 | Independent integration of the angular equilibrium ODE checks the closed form and face traction. This is a sound check despite the separate reflection error. |
| Verification shows numerical results | 5 / 5 | Rerun reproduces max angular discrepancy 1.554 × 10⁻¹⁵ and g(π) = 1.595 × 10⁻¹⁶. The n = 1 limiting calculation also reproduces I_J = 0.7853981633972231 versus π/4. |
| Verification interprets agreement | 5 / 5 | Explains that g(π) should vanish for the traction-free face and that the n = 1 results agree with the known linear solution. |
| Missing fracture physics | 5 / 5 | Correctly says that initiation requires a separate failure criterion, such as critical J or critical slip at a characteristic distance. |

**Total: 99 / 100.** Category totals: model 29/30; numerical correctness 40/40; completion and conventions 10/10; verification 15/15; missing fracture physics 5/5.

The submission scripts and reference were rerun from copies in `scratch/`. The reference's two independent routes and its assertions pass; no reference error was found, so an alternative grade is unnecessary. Detailed audit results are in [scratch/audit_output.txt](<matched-runs>/grading/M43/scratch/audit_output.txt).

Scientific errors:

- The stress reflection in `submission/solve.py:79`–`83` and `:181`–`185` is inconsistent with odd displacement. Correct polar parity is odd τ_rz and even τ_θz; the code uses even τ_rz and odd τ_θz. At r = 0.0200 mm, θ = −60°, `tau_cart` gives Cartesian stresses (−100.2274, −102.4223) MPa, whereas its own odd displacement and the constitutive law require (+100.2274, +102.4223) MPa. At θ = 0, `np.sign(0)` also incorrectly sets τ_yz to zero instead of approximately +153.1222 MPa. This single implementation error accounts for the symmetry deduction; the correctly evaluated upper-half outputs retain full numerical credit.
- The opening code docstring has stale formulas: line 9 puts a minus sign in the J radial exponent, which should be 1 + (α_w−1)(1+1/n); lines 15–16 give angular functions incompatible with the corrected expressions actually used by the executable code. The answer's finite-flux argument and operative angular formulas are correct.

Unsupported claims:

- The claimed verification of the “full field” is broader than the checks performed. Cartesian equilibrium samples and the displacement-gradient comparison all lie in the upper half-plane, so they do not detect the incompatible lower-half stresses or the jump across the positive x-axis. The independent angular ODE check still satisfies all three verification criteria.
- The stated 1.0–1.5% engineering-shear range at requested positions is imprecise: the requested 60° stress point has |γ| = 0.01538553, or 1.538553%. At the face-slip positions the shear is below 1%. This does not change the requested calculations.

Presentation problems:

- The symbol p denotes the p-Laplacian exponent 9/7 and then the stress exponent 2/9 without explicitly distinguishing them.
- “θ(π)=0-condition” should identify zero face traction, τ_θz(π)=0.
- The ODE verification is called “stiff IVP integration,” but the supplied `solve_ivp` call uses its default explicit RK45 method; the numerical check itself succeeds.

Tolerance: **All seven requested numbers, including both signed stresses, are within the stated 2% relative tolerance; the largest error in the answer table is 0.015465%.**
