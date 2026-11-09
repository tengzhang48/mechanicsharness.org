| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: leading tensor compliance | 10 / 10 | Correctly retains the elastic volume response, derives effective modulus T and ν_eff = 1/2 − (1/2 − ν)T/E, and permits nonzero ε_zz under plane stress. |
| Model: traction-free tensile leading field | 10 / 10 | Correct mode-I Williams field after the constitutive reduction, with r^−1/2 stresses, r^1/2 displacements, and the correct polar angular functions. The yield-intercept strain is correctly identified as bounded. |
| Model: J normalization | 5 / 5 | Uses the appropriate plane-stress local relation J = K_tip²/T and independently verifies its flux. |
| Model: transverse inversion | 5 / 5 | Relates ε_zz to both in-plane stresses, obtains C_zz = −P√(JT/(2π)), and selects the unique admissible root 0 < T < E. |
| Numerical: T | 12 / 12 | Reports 4816 MPa, with 4816.47 MPa also supplied, versus 4816.472187 MPa. Relative errors are 0.009804% and 0.00004541%, respectively; both satisfy 1%. |
| Numerical: C_rr at 60° | 5 / 5 | +0.003942 mm^1/2 versus +0.003942048238; relative error 0.001224%, below 2%. |
| Numerical: C_θθ at 60° | 5 / 5 | +0.0006332 mm^1/2 versus +0.0006331781351; relative error 0.003453%, below 2%. |
| Numerical: C_γrθ at 60° | 5 / 5 | +0.005731 mm^1/2 versus +0.005731131133; relative error 0.002288%, below 2%. Correct engineering shear and positive sign. |
| Numerical: C_y at 90° | 5 / 5 | +0.009139 mm^1/2 versus +0.009139018941; relative error 0.0002073%, below 2%. |
| Numerical: lower radius endpoint | 4 / 4 | Explicitly reports 0.0100 mm for the nonempty leading-field interval. Examines all three principal strains, including the controlling ε_zz; agrees with the reference to rounding. |
| Numerical: upper radius endpoint | 4 / 4 | Explicitly reports 0.0613 mm versus 0.06132522855 mm; relative error 0.04114%, below 2%. Correctly intersects the stress bound with r ≤ 0.200 mm. |
| Completion: seven outputs | 5 / 5 | T, the three angular strain coefficients, displacement coefficient, and both endpoints are all identifiable. The final written interval is explicitly nonempty and includes both endpoints. |
| Conventions: units | 2 / 2 | Correct MPa, mm^1/2, and mm units. |
| Conventions: tensor/engineering shear | 2 / 2 | Explicitly defines C_γrθ using 2ε_rθ and includes the factor of two in the calculation. |
| Conventions: polar/displacement signs | 1 / 1 | Polar rotation, positive shear, and positive u_y at 90° agree with the specified basis and mode-I opening. |
| Verification: independent condition | 5 / 5 | solve.py directly integrates W N_x − t_i u_i,x using numerical displacement derivatives; verify2.py independently differentiates and rotates the displacement strain. These checks go beyond recalibration of C_zz. |
| Verification: quantitative evaluation | 5 / 5 | Reruns give J = 0.80000000 N/mm on all three tested circles. Finite differences give (0.003942048, 0.000633178, 0.005731131), matching the analytic strain coefficients to the displayed nine decimal places. |
| Verification: interpretation | 5 / 5 | Correctly explains that flux agreement validates J–K–T normalization and displacement/strain agreement validates the angular field. These reproducible checks suffice for all 15 verification points. |
| Scope: meaning of screens | 3 / 3 | Explicitly identifies the screens as checks on the leading field, not finite-r error estimates or failure criteria; acknowledges limits of physical applicability and distinguishes the unreduced-law comparison. |
| Scope: crack initiation | 2 / 2 | States that initiation needs an additional toughness, fracture, or damage criterion absent from the supplied data. |

**Total: 100 / 100.** Category totals: model 30/30; numerical correctness 40/40; completion/conventions 10/10; verification 15/15; scope 5/5.

Scientific errors:

- No errors in the scored leading-field derivation or requested numbers; the reference rerun also passed its tensor, flux, compatibility, and limiting checks, with no reference error identified.
- The ancillary claim that this local leading field necessarily presumes small-scale yielding is too restrictive. The specified deformation law admits the stated high-stress local reduction without that additional global loading assumption. This does not change the correctly derived field or trigger a rubric deduction.

Unsupported claims:

- The proposed blunting estimate δ_t ≈ d_n J/σ_y with d_n ≈ 0.3–0.8 is not derived for the supplied material or established by the problem. Its comparison with the inner screen is speculative, as the submission itself acknowledges.
- Calling the reduced law “accurate” solely because σ_e ≥ 10σ_y is stronger than the screen establishes: the discarded factor σ_y/σ_e can still be 0.1. The answer nevertheless explicitly denies a finite-r error guarantee, meeting the scope criterion.

Presentation problems:

- solve.py computes the correct leading lower bound, but assigns its final printed interval to r_b_full = 0.00926930797 mm instead. ANSWER.md clearly selects the correct leading interval [0.0100, 0.0613] mm and labels the alternative calculation, so the final-answer endpoints receive full credit.
- The extra blunting and applicability speculation adds unnecessary ambiguity to an otherwise clear answer.

Tolerance: all seven requested final-answer numbers pass their stated relative tolerances (T: 1%; the other six: 2%); the largest reported relative error is 0.04114%.
