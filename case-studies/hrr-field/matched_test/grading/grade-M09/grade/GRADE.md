# Grade

| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: leading tensor compliance and transverse strain | 10 / 10 | Correct coefficients α1 and α2 retain the elastic volume response: α1 − 3α2 = (1 − 2ν)/E. The equivalent modulus is T and ν2 = 1/2 − (1/2 − ν)T/E; ε_zz = −α2(σ_xx + σ_yy) remains free and nonzero. The incorrect descriptive word “anisotropic” does not affect the actual reduction. |
| Model: traction-free tensile leading field | 10 / 10 | Correct classical mode-I angular stresses, r^(−1/2) stress/strain behavior, and r^(1/2) displacement behavior. Justified by the homogeneous isotropic linear leading compliance and traction-free faces. |
| Model: J normalization | 5 / 5 | Uses K = √(JT), hence J = K²/T, and independently evaluates the specified contour flux. |
| Model: transverse inverse calibration | 5 / 5 | Correct equation C_zz = −2α2√(JT/(2π)); the bracketed root solve gives the unique admissible T = 4816.472187302739 MPa. The calibration function is strictly increasing on 0 < T < E. |
| Numerical: T | 12 / 12 | Reports 4816.472 MPa versus 4816.472187302739 MPa: relative error 0.000003889%, below 1%. The coarser headline value 4816 MPa also passes, with error 0.009804%. |
| Numerical: C_rr(60°) | 5 / 5 | Reports +0.003942 mm^1/2 versus +0.003942048237733179: relative error 0.001224%, below 2%. |
| Numerical: C_θθ(60°) | 5 / 5 | Reports +0.0006332 mm^1/2 versus +0.0006331781350710297: relative error 0.003453%, below 2%. |
| Numerical: C_γrθ(60°) | 5 / 5 | Reports +0.005731 mm^1/2 versus +0.005731131133456488: relative error 0.002288%, below 2%; correctly includes the factor of two for engineering shear. |
| Numerical: C_y(90°) | 5 / 5 | Reports +0.009139 mm^1/2 versus +0.009139018941480233: relative error 0.0002073%, below 2%. Both displacement and radial-strain routes reproduce the result. |
| Numerical: lower radius endpoint | 4 / 4 | Reports 0.01000 mm versus 0.01000000000 mm: zero error to displayed reference precision. Examines all three principal strains, including the controlling ε_zz, and explicitly identifies a nonempty interval. |
| Numerical: upper radius endpoint | 4 / 4 | Reports 0.06133 mm versus 0.06132522854990914 mm: relative error 0.007781%, below 2%. Correctly intersects the effective-stress bound with r ≤ 0.200 mm; the interval is explicitly nonempty. |
| Completion: seven identifiable outputs | 5 / 5 | T, three strain coefficients, the displacement coefficient, and both radius endpoints are all reported. |
| Conventions: units, shear, and signs | 5 / 5 | Correct MPa, mm, and mm^1/2 units (2); correct tensor/engineering shear distinction (2); consistent polar basis and positive requested displacement/shear signs (1). |
| Verification: independent condition | 5 / 5 | The submitted solver computes J directly from W = σ:ε/2 and Cartesian displacement derivatives on circular contours. This checks the field normalization independently of substitution into C_zz. |
| Verification: quantitative evaluation | 5 / 5 | Rerun gives J = 0.800000000000475, 0.799999999994723, 0.799999999999604, and 0.799999999985885 N/mm at R = 0.1, 0.5, 1.0, and 2.0 mm. Maximum absolute discrepancy from 0.800 is 1.42 × 10^(−11) N/mm. |
| Verification: interpretation | 5 / 5 | Correctly interprets the contour agreement as checking J = K²/T. These are calculations of the formal leading field; they do not extend the stated finite-radius matching range. This one reproduced check suffices for all 15 verification points. |
| Scope: meaning of the screening interval | 3 / 3 | Describes the screens as consistency checks rather than error bars, acknowledges the approximately 10% omitted constitutive factor at the high-stress boundary, and uses the requested leading field. Does not equate the numerical tolerance with finite-radius physical accuracy. |
| Scope: initiation and missing physics | 2 / 2 | Explains that J and the constitutive parameters do not determine initiation without a fracture/failure/damage criterion or toughness. |

**Total: 100 / 100.**

The copied submission solver reproduces the reference results. The reference solver also reruns successfully with all its assertions passing; no reference error was found. Detailed rerun outputs and tolerance comparisons are saved in `scratch/`.

(a) Scientific errors:

- The phrase “anisotropic-in-shear linear elastic law” is incorrect: the leading tensor law is isotropic, with different effective shear and bulk responses. The subsequent equations and effective-modulus identification are correct, so this wording slip does not lose rubric credit.
- The asymptotic explanation incorrectly says “1 − σ_y/σ_e → 0.” That factor tends to 1; the omitted relative correction σ_y/σ_e tends to 0, and its strain contribution is bounded. The earlier limiting argument and implemented leading law use the correct limit, so this prose error does not affect the scored calculation.

(b) Unsupported claims:

- The additional reported equilibrium, symmetry, and compatibility PASS residuals are not reproduced by the self-contained submitted files: the auxiliary verifier fails to import its checker. They receive no separate evidentiary credit. The independently reproduced direct-J check already satisfies the full verification rubric.

(c) Presentation problems:

- No material presentation problem: all requested quantities, units, signs, and the nonempty closed interval are identifiable. The two misleading phrases are recorded above.

Tolerance verdict: all seven requested numbers pass—T is within 1% relative tolerance, and each of the other six numbers is within 2%.
