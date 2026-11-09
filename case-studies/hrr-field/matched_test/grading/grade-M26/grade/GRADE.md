# Grade

| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: leading tensor compliance | 10 / 10 | Correctly removes the bounded yield-intercept contribution, retains the elastic bulk modulus, derives E_eq = T and ν_eq = 1/2 − (1/2 − ν)T/E, and permits ε_zz = −ν_eq(σ_xx + σ_yy)/T. |
| Model: tensile leading field | 10 / 10 | Justifies the traction-free, symmetric mode-I Williams field after reduction, with stresses proportional to r^(−1/2), displacements proportional to √r, and correct angular functions. |
| Model: J normalization | 5 / 5 | Uses the correct plane-stress relation J = K_I²/T; the directly integrated flux tends to 0.800 N/mm. |
| Model: transverse inverse | 5 / 5 | Uses C_zz = −2ν_eq√(J/(2πT)) and solves for admissible T = 4816.47 MPa. This equation is monotonic on 0 < T < E and has a unique admissible root. |
| Numerical: T | 12 / 12 | 4816.47 MPa; relative error 0.00004541%, within 1%. |
| Numerical: C_rr(60°) | 5 / 5 | +0.0039420 mm^1/2; relative error 0.001224%, within 2%. |
| Numerical: C_θθ(60°) | 5 / 5 | +0.00063318 mm^1/2; relative error 0.0002945%, within 2%. |
| Numerical: C_γrθ(60°) | 5 / 5 | +0.0057311 mm^1/2; relative error 0.0005432%, within 2%; includes the engineering-shear factor of two. |
| Numerical: C_y(90°) | 5 / 5 | +0.0091390 mm^1/2; relative error 0.0002073%, within 2%. |
| Numerical: lower radius endpoint | 4 / 4 | 0.010000 mm; agrees to roundoff. Examines all three principal strains, with the out-of-plane coefficient −0.005 binding. Explicitly reports a nonempty interval. |
| Numerical: upper radius endpoint | 4 / 4 | 0.061325 mm; relative error 0.0003727%, within 2%. The coarser 0.0613 mm also passes (0.04114% error). Correctly combines the effective-stress bound with r ≤ 0.200 mm. |
| Completion: seven outputs | 5 / 5 | All seven requested numbers are reported and identifiable, including both interval endpoints. |
| Conventions: units | 2 / 2 | All requested outputs have the correct MPa, mm^1/2, or mm units. An auxiliary unit typo is noted below. |
| Conventions: tensor/engineering shear | 2 / 2 | Explicitly computes C_γrθ for 2ε_rθ, while using tensor shear in the strain law. |
| Conventions: polar/displacement signs | 1 / 1 | Polar rotation agrees with the stated basis; the requested coefficients and u_y have the correct signs, with no additive translation. |
| Verification: independent condition | 5 / 5 | Integrates the polar strains and checks the independent shear compatibility equation U_r′ − U_θ/2 = 2C_εrθ. This is not a substitution into C_zz. |
| Verification: quantitative evaluation | 5 / 5 | The unmodified core reproduces a 7.78×10^−13 compatibility residual and C_y = 0.009139019 by strain integration, closed displacement formula, and the vertical line integral. |
| Verification: interpretation | 5 / 5 | Correctly interprets those agreeing displacement calculations as verification. This complete check earns the 15 points; the erroneous optional Cartesian shear diagnostic is not credited. |
| Scope: meaning of screens | 3 / 3 | Identifies the high-stress, three-principal-strain, and matching screens as defining a working interval and expressly says they are not finite-r error estimates. |
| Scope: crack initiation | 2 / 2 | States that initiation requires an additional toughness, critical-condition, or damage criterion absent from the supplied law. |

**Total: 100 / 100.** No error was found in the reference, so no alternative grading is required. The issues below do not invalidate any requested result or the separate, successful verification that satisfies the rubric.

Scientific errors:

- `submission/solve.py:129` labels its diagnostic as 1/2(∂u_x/∂y + ∂u_y/∂x), but actually uses ∂u_y/∂y in the second term. Its log therefore shows shear residuals 0.0132, 0.0182, and 0.0514. This is an error in the diagnostic, not in the reported shear coefficient or displacement field.

Unsupported claims:

- The assertion in `ANSWER.md:72` that the submitted Cartesian finite-difference strain check agrees to approximately 10^−12 is contradicted by its own shear diagnostic and log. The separate polar compatibility check does support the awarded verification credit. Correcting the derivative in a grader-only audit confirms the source of the discrepancy; that correction is not counted as the solver's verification.
- The finite-circle calculation combines full-law W with leading-field stresses and displacements. Its reproduced J values verify the limiting normalization, but do not measure corrections or error in the complete finite-r solution. Thus the quoted <0.06% flux discrepancy cannot establish such an error bound.
- The additional claim that thickness must be known to validate the imposed plane-stress model introduces a condition outside the mathematical problem, which explicitly prescribes that constraint. Matching corrections being negligible also does not quantify every correction from the constitutive yield intercept.

Presentation problems:

- The auxiliary quantity K_I/√(2π) is labeled MPa in both the answer and log; its unit is MPa·mm^1/2. Units of the seven requested outputs are correct.
- The verification prose does not acknowledge the large shear residuals in the accompanying log.

Tolerance verdict: **All seven requested numbers are within tolerance: T within 1%, and the other six within 2%; the screened interval is explicitly nonempty.**
