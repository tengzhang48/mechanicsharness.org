# Grade

Graded against RUBRIC.md and the supplied reference. Unmodified copies of both scripts ran successfully in `scratch/`. The submission's computed outputs agree with the reference to a maximum relative difference of 5.55 × 10⁻¹⁶. The numerical comparisons below use the rounded results in the answer's opening table. No error was found in the reference.

| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: high-stress tensor compliance | 10 / 10 | Correctly obtains effective modulus T and ν* = 1/2 − (1/2 − ν)T/E, preserves the elastic bulk modulus, and allows nonzero ε_zz. Original and reduced bulk moduli both equal 58,333.3333 MPa. The erroneous finite-stress expression noted below does not affect this leading reduction. |
| Model: traction-free tensile leading field | 10 / 10 | Justifies the plane-stress mode-I Williams field after the constitutive reduction; the implemented angular stresses, r⁻¹ᐟ² stresses/strains, and √r displacements are correct. |
| Model: normalization using J | 5 / 5 | Uses J = K²/T and explains why bounded strain corrections vanish from the shrinking-contour leading flux. The contour calculation verifies this normalization. |
| Model: transverse inverse calibration | 5 / 5 | Correctly relates ε_zz to the in-plane stress trace and solves the monotone inverse equation for an admissible unique T. The equation immediately before the displayed typo and the implemented equation are correct. |
| Numerical: T | 12 / 12 | 4816 MPa; relative error 0.009804%, within 1%. Code gives 4816.472187 MPa. |
| Numerical: C_rr(60°) | 5 / 5 | +0.003942 mm^1/2; relative error 0.001224%, within 2%. |
| Numerical: C_θθ(60°) | 5 / 5 | +0.0006332 mm^1/2; relative error 0.003453%, within 2%. |
| Numerical: C_γrθ(60°) | 5 / 5 | +0.005731 mm^1/2; relative error 0.002288%, within 2%. Correct engineering shear factor and positive sign. |
| Numerical: C_y(90°) | 5 / 5 | +0.009139 mm^1/2; relative error 0.0002073%, within 2%. |
| Numerical: lower radius endpoint | 4 / 4 | 0.0100 mm; relative error effectively zero. Examines all three principal strains: at this endpoint they are (+0.0264151, +0.0264151, −0.0500000), so ε_zz controls. Explicitly recognizes a nonempty interval. |
| Numerical: upper radius endpoint | 4 / 4 | 0.0613 mm; relative error 0.041139%, within 2%. Correct effective-stress bound, combined with r ≤ 0.200 mm; explicitly recognizes a nonempty interval. |
| Completion: seven identifiable outputs | 5 / 5 | Reports T, four coefficients, and both radius endpoints. |
| Conventions: units | 2 / 2 | All seven requested outputs carry correct units. The auxiliary stress-coefficient unit typo is noted below. |
| Conventions: tensor/engineering shear | 2 / 2 | Uses tensor ε_xy in the compliance and doubles ε_rθ for the requested engineering shear. |
| Conventions: polar/displacement signs | 1 / 1 | Correct polar rotation and tensile displacement signs; no additive rigid translation. |
| Verification: independent condition | 5 / 5 | Implements the actual J contour integral using energy, tractions, and finite-difference displacement gradients, independently of substituting into C_zz. |
| Verification: quantitative evaluation | 5 / 5 | Rerun gives J = 0.8000000000082, 0.7999999999963, and 0.8000000000094 N/mm at r = 0.001, 0.010, and 0.100 mm. |
| Verification: interpretation | 5 / 5 | Correctly interprets agreement with supplied J as validation of the leading-field normalization; the script also provides displacement/strain compatibility checks. One complete contour check earns all 15 verification points. |
| Scope: meaning of screens | 3 / 3 | Applies the screens to the stated leading field and explicitly distinguishes them from full finite-r error estimates. Does not equate the numerical tolerance with physical accuracy. |
| Scope: crack initiation | 2 / 2 | Explains that initiation needs an additional fracture, failure, or damage criterion absent from the supplied data. |

**Total: 100 / 100.**

The following issues concern ancillary explanations or transcription; the answer and code still satisfy every scored requirement.

## (a) Scientific errors

- In `submission/ANSWER.md`, line 24, the stress potential that differentiates to strain is complementary energy Φ(σ), rather than the problem's strain work W = ∫σ:dε. Above yield, Φ = Φ_el + (A/2)(σ_e − σ_y)², whereas W = Φ_el + (A/2)(σ_e² − σ_y²), with A = 1/T − 1/E and Φ_el = [(1+ν)σ:σ − ν(tr σ)²]/(2E). For uniaxial stress 30 MPa, these are 0.0450956 and 0.0837626 MPa, respectively. They share the same leading quadratic term, and the submitted leading-field J calculation uses the correct W = σ:ε/2.
- The finite-stress ε_zz expression on that same line puts the yield multiplier outside the entire elastic-plus-plastic coefficient. The correct expression is ε_zz = −[ν/E + (A/2)⟨1 − σ_y/σ_e⟩](σ_xx + σ_yy): only the plastic contribution receives that multiplier. The stated leading limit and the finite-stress check in the code are correct.

## (b) Unsupported claims

- The claimed separately executed polar-derivative check at θ = 60° is not documented by the supplied script. It prints stress-derived polar coefficients at 60° and checks Cartesian displacement gradients at four other angles. This does not reduce verification credit because its fully documented contour check suffices.
- The uncertainty statement supplies no input uncertainty bounds or propagation calculation. Printed significant figures alone do not establish an uncertainty estimate for T. This does not affect the specified numerical-tolerance grading.

## (c) Presentation problems

- The displayed calibration equation on line 30 has √T in the denominator where T is required: C_zz = −2ν*√(JT)/(T√(2π)). Its literal printed form evaluates to −0.347004 rather than −0.005000 and has the wrong dimensions. The preceding equation and code correctly implement the intended calibration.
- The auxiliary K/√(2π) on line 14 should have units MPa·mm^1/2, rather than MPa. Requested-output units are correct.
- The angular stress amplitude C in the prose is not explicitly defined, and u denotes both displacement and θ/2 in the same paragraph.

Tolerance verdict: all seven requested numbers are within tolerance—T within 1%, and each of the other six within 2%; the reported radius interval is explicitly nonempty and includes both endpoints.

Reproducible evidence: `scratch/submission_run.txt`, `scratch/reference_run.json`, and `scratch/grader_checks.json` (generated by `scratch/grader_checks.py`).
