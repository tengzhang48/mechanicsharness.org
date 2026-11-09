| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: high-stress tensor compliance | 10 / 10 | Correctly obtains the isotropic pair T and ν_t = 1/2 − (1/2−ν)T/E. The deviatoric plastic term leaves the elastic volume response intact; ε_zz is free and nonzero. |
| Model: tensile traction-free leading field | 10 / 10 | Correct plane-stress mode-I field, justified by the reduced linear operator and symmetry, with r^−1/2 stresses and r^1/2 displacements. The final reconstructed displacement is compatible with that field. |
| Model: J normalization | 5 / 5 | Uses J = K_t²/T and K_t = √(JT); the direct contour calculation in the final verify3.py validates the factor. |
| Model: transverse inverse measurement | 5 / 5 | Correct equation C_zz(0) = −2ν_t√[J/(2πT)] and admissible unique root T = 4816.472187 MPa. The prose reverses the direction of monotonicity, but the equation and inversion are correct. |
| Numerical: T | 12 / 12 | Reports 4816.5 MPa versus 4816.472187302739 MPa; relative error 0.00057745%, below 1%. |
| Numerical: C_rr(60°) | 5 / 5 | Reports +0.0039420 versus +0.003942048237733179 mm^1/2; relative error 0.00122367%, below 2%. |
| Numerical: C_θθ(60°) | 5 / 5 | Reports +0.00063318 versus +0.0006331781350710297 mm^1/2; relative error 0.00029453%, below 2%. |
| Numerical: C_γrθ(60°) | 5 / 5 | Reports +0.0057311 versus +0.005731131133456488 mm^1/2; relative error 0.00054323%, below 2%. Correct engineering shear factor and sign. |
| Numerical: C_y(90°) | 5 / 5 | Reports +0.0091390 versus +0.009139018941480233 mm^1/2; relative error 0.00020726%, below 2%. The final displacement calculation reproduces this value. |
| Numerical: lower radius endpoint | 4 / 4 | Reports 0.01000 mm, agreeing with 0.01000000000 mm to rounding. Examines all three principal strains and correctly finds that ε_zz controls. Explicitly identifies the interval as nonempty. |
| Numerical: upper radius endpoint | 4 / 4 | Reports 0.06133 versus 0.06132522854990914 mm; relative error 0.00778057%, below 2%. Correct effective-stress screen; the 0.200 mm bound is nonbinding. Interval is nonempty. |
| Completion: seven outputs | 5 / 5 | T, the three angular strain coefficients, the displacement coefficient, and both radius endpoints are reported and identifiable. |
| Conventions: units, shear, signs | 5 / 5 | Correct MPa/mm^1/2/mm units (2), tensor versus engineering shear convention (2), and polar/displacement signs (1). |
| Verification: independent condition | 5 / 5 | verify3.py directly evaluates the J flux from work density, traction, and displacement derivatives, independently checking the normalization. It also checks sym(∇u) against the leading tensor law. |
| Verification: quantitative evaluation | 5 / 5 | Rerun gives leading-field J = 0.7999999999916344 N/mm at r = 0.01 mm. Re-evaluation of the submitted strain checks at six angles gives maximum absolute mismatch 2.85×10^−11. |
| Verification: interpretation | 5 / 5 | Correctly connects flux agreement to J = K_t²/T and strain agreement to the angular coefficients and compatible displacement. This valid final check earns full credit despite the ancillary errors listed below. |
| Scope: meaning of screens | 3 / 3 | Explicitly treats the screens as checks of the leading field, distinguishes them from full-solution error estimates, and explains the small-strain and matching restrictions. |
| Scope: crack initiation | 2 / 2 | Correctly states that the damage-free constitutive/loading data need an additional fracture or failure criterion to determine initiation. |

**Total: 100 / 100.** The final leading-field results meet every scored criterion. The rubric permits one valid independent verification to earn all 15 verification points; the unsuccessful earlier checks do not negate the valid final check. The reference rerun passes its assertions and reproduces the recorded key; no reference correction or alternative score is warranted.

**(a) Scientific errors**

- ANSWER.md lines 37 and 60 describe the signed calibration coefficient as strictly decreasing. With b = (1−2ν)/E, C_zz(T) = −√[J/(2π)](T^−1/2 − bT^1/2), whose derivative is positive. Its positive magnitude decreases. The uniqueness conclusion and computed root remain correct.
- ANSWER.md line 53 incorrectly rejects the recalled displacement formula. Both of its brackets equal κ−cosθ, and with the effective moduli it is the same field as the final reconstruction. The two implementations agree to 4.34×10^−19 mm on the sampled angles. The apparent J = 0.0885321012 N/mm failure comes from solve.py lines 123–126 using u_y,y instead of u_y,x. Correcting only that derivative in the scratch check gives 0.7999999999962559 N/mm.
- The retained verify2.py reconstruction is incorrect: its rerun gives C_y = −0.000555901854 mm^1/2, J = −0.09549925 N/mm, and large strain discrepancies. The final answer instead uses the correct verify3.py reconstruction, so these earlier results are not the scoring target.
- The complete-law work formula displayed in ANSWER.md line 57 applies to yielded points. Below σ_e = σ_y the plastic work is zero; the factor needs the yield cutoff. The submitted code correctly includes that cutoff.

**(b) Unsupported claims**

- The asserted contour-error rate “∝ √ρ” in ANSWER.md line 57 is contradicted by the rerun: the full-work-density correction is 1.64538×10^−4 N/mm at ρ = 0.1 mm and 5.20460×10^−6 N/mm at ρ = 0.01 mm. Their ratio is 31.614, indicating approximately ρ^1.5 for this calculation. Convergence and the leading-field normalization are nevertheless verified.
- The statement that the recalled displacement formula passes strain checks only on shear-free rays is unsupported and contradicted by its equivalence to the compatible final field.

**(c) Presentation problems**

- solve.py and verify2.py retain conflicting failed checks alongside the successful verify3.py; running every supplied script produces inconsistent verification outputs. The final answer identifies verify3.py, but the unsuccessful code and mistaken explanation of the earlier failure make the evidence harder to follow.
- The abandoned-calculation narrative is unnecessary for evaluating the final mechanics. The compact Cartesian stress notation also compresses three components into two expressions with a paired sign, reducing readability.

All seven requested numbers are within tolerance: T within 1%, and each of the other six numbers within 2%.
