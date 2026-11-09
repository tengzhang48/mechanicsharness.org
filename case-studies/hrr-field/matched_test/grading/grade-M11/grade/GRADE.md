# Grade

Graded against `RUBRIC.md` and the recorded reference. Both scripts were copied into `scratch/` and rerun successfully. The reference assertions passed; no reference discrepancy was found. Numerical errors below conservatively use the coarsest rounding in the submission's main results table. Its code reproduces all seven reference targets to better than 10⁻¹² relative error.

| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: high-stress tensor compliance | 10 / 10 | Correctly retains the elastic bulk modulus, derives the modified shear compliance and equivalent modulus T, and allows nonzero ε_zz. The resulting ν_eff = 0.486238650893 agrees with ν_t = 1/2 − (1/2 − ν)T/E. |
| Model: traction-free tensile leading field | 10 / 10 | Uses the appropriate Mode-I field with r⁻¹ᐟ² stresses and √r displacements. Explains that the yield offset contributes bounded strains and gives the correct angular stresses and displacement field. |
| Model: normalization using J | 5 / 5 | Uses the plane-stress relation J = K_I²/E_eff = K_I²/T, with a numerical contour-flux check. |
| Model: inverse transverse calibration | 5 / 5 | Relates ε_zz to the sum of the in-plane normal stresses and solves the correct inverse equation. The scalar inversion gives the admissible solution T = 4816.472187 MPa, with 0 < T < E. |
| Numerical: T | 12 / 12 | Reported 4816 MPa (also 4816.47); relative error 0.009804%, within 1%. |
| Numerical: C_rr at 60° | 5 / 5 | Reported +0.003942 mm^1/2; relative error 0.001224%, within 2%. |
| Numerical: C_θθ at 60° | 5 / 5 | Reported +0.0006332 mm^1/2; relative error 0.003453%, within 2%. |
| Numerical: C_γrθ at 60° | 5 / 5 | Reported +0.005731 mm^1/2; relative error 0.002288%, within 2%. Correct engineering shear and sign. |
| Numerical: C_y at 90° | 5 / 5 | Reported +0.009139 mm^1/2; relative error 0.0002073%, within 2%. |
| Numerical: lower radius endpoint | 4 / 4 | Reported 0.0100 mm; relative error effectively zero, within 2%. Examines all three principal strains, including the controlling ε_zz, and explicitly identifies the interval as nonempty. |
| Numerical: upper radius endpoint | 4 / 4 | Reported 0.0613 mm; relative error 0.04114%, within 2%. Correct stress-invariant screen and intersection with r ≤ 0.200 mm; explicitly nonempty. |
| Completion: seven identifiable outputs | 5 / 5 | Reports T, all four requested strain/displacement coefficients, and both radius endpoints. |
| Conventions: units | 2 / 2 | Correct MPa, mm^1/2, and mm units for the requested outputs. |
| Conventions: tensor and engineering shear | 2 / 2 | Consistently distinguishes ε_rθ from 2ε_rθ, including in the rotation, constitutive calculation, and kinematic check. |
| Conventions: polar and displacement signs | 1 / 1 | Rotation and displacement formulas use the stated basis and give the correct positive coefficients. |
| Verification: independent condition | 5 / 5 | Actually computes strains from differentiated displacements and evaluates J from W and u_,x. These test compatibility and normalization independently of substitution into C_zz. |
| Verification: quantitative evaluation | 5 / 5 | Reproduced maximum relative strain-coefficient discrepancy 7.56 × 10⁻⁷, confirming the claimed <10⁻⁶. Numerical J values at r = 0.001, 0.01, 0.1 mm are 0.800000000033, 0.800000000059, 0.799999999997 N/mm. |
| Verification: interpretation | 5 / 5 | Correctly interprets strain agreement as a kinematic/constitutive check and the contour results as confirmation of J normalization and radius independence for the reconstructed leading field. |
| Scope: meaning of the screens | 3 / 3 | Explicitly calls them defined checks of the leading field rather than finite-r error estimates, explains each bound, and identifies the working interval. Does not equate numerical tolerance with physical accuracy. |
| Scope: missing initiation physics | 2 / 2 | Explains that an additional fracture, damage, cohesive, or critical-opening criterion is needed; a supplied J alone does not determine initiation. |

**Total: 100 / 100** (model 30/30; numerical 40/40; completion/conventions 10/10; verification 15/15; scope 5/5).

**(a) Scientific errors**

- None found in the requested leading-field calculation, constitutive reduction, screening interval, or initiation conclusion.

**(b) Unsupported claims**

- None material to the rubric. The submitted code reproduces the claimed independent kinematic and contour-integral checks; credit does not rely on the calibration substitution alone.

**(c) Presentation problems**

- Minor wording in §6: saying measurement errors “propagate ~quadratically” conflates inverse-square parameter dependence with error propagation. At fixed J, small fractional changes obey δT/T ≈ −1.89 δ|C_zz|/|C_zz| at this solution. This optional sensitivity wording does not change any rubric score.

**Tolerance verdict:** All seven requested numbers are within tolerance: T within 1%, and the other six within 2%; the reported interval is nonempty.
