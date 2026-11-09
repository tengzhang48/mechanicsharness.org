| rubric item | points awarded / max | reason |
|---|---:|---|
| Leading constitutive compliance | 10 / 10 | Correctly collects A and B, obtains effective Young modulus T and ν* = 1/2 − (1/2 − ν)T/E, and permits nonzero ε_zz under plane stress. The deviatoric plastic response preserves the elastic volume response. |
| Tensile leading field and justification | 10 / 10 | Correct reduced isotropic mode-I plane-stress field, traction-free angular functions, r^−1/2 stresses, and compatible r^1/2 displacements. |
| J normalization | 5 / 5 | Correctly uses J = K²/T and K = √(JT). This analytical normalization is correct despite the flawed supplementary numerical J routine. |
| Transverse inverse and admissible T | 5 / 5 | Correctly relates ε_zz to the in-plane stress trace and solves C_zz(0) = −2ν*√(J/(2πT)). The supplied code brackets the admissible root in 0 < T < E; this inverse relation has a unique root there. |
| T | 12 / 12 | Reports 4816.47 MPa versus 4816.472187302739 MPa; relative error 0.000045413%, below 1%. |
| C_rr at 60° | 5 / 5 | Reports +0.0039421 versus +0.003942048237733179 mm^1/2; relative error 0.00131308%, below 2%. |
| C_θθ at 60° | 5 / 5 | Reports +0.00063318 versus +0.0006331781350710297 mm^1/2; relative error 0.000294535%, below 2%. |
| C_γrθ at 60° | 5 / 5 | Reports +0.0057311 versus +0.005731131133456488 mm^1/2; relative error 0.000543234%, below 2%. Correct engineering shear factor and sign. |
| C_y at 90° | 5 / 5 | Reports +0.0091390 versus +0.009139018941480233 mm^1/2; relative error 0.000207259%, below 2%. |
| Lower radius endpoint | 4 / 4 | Reports 0.01000 mm, agreeing with the reference. Examines all three principal strains and correctly identifies ε_zz as controlling. Explicitly reports a nonempty interval. |
| Upper radius endpoint | 4 / 4 | Reports 0.06133 versus 0.06132522854990914 mm; relative error 0.00778057%, below 2%. Correct effective-stress screen and inactive 0.200 mm cap; interval explicitly nonempty. |
| All requested outputs | 5 / 5 | All seven requested numerical outputs are identifiable, including both interval endpoints. |
| Units | 2 / 2 | Correct MPa, mm^1/2, and mm units. |
| Tensor/engineering shear convention | 2 / 2 | Correctly distinguishes ε_rθ from 2ε_rθ in formulas and code. |
| Polar/displacement signs | 1 / 1 | Consistent polar rotation, positive requested shear coefficient, and positive upper-half-plane displacement. |
| Verification: independent condition | 5 / 5 | The finite-difference sym(∇u) comparison with stress-law strains is a genuine independent compatibility check, separate from fitting C_zz. No credit is attributed to the erroneous direct J check. |
| Verification: quantitative evaluation | 5 / 5 | Unchanged solve.py reproduces the compatibility comparison at r = 0.1 mm and θ = 10°, 60°, 90°, 150°. Maximum absolute strain difference is 5.60591 × 10^−12, supporting the stated agreement to approximately 10^−10. |
| Verification: interpretation | 5 / 5 | Correctly interprets the successful strain/displacement comparison as confirming consistency of the reduced mode-I field. The rubric permits one successful independent check to earn all 15 verification points. |
| Scope of screening approximation | 3 / 3 | Applies the screens to the leading field, explicitly denies finite-r error bounds, and identifies the resulting working interval. Some validity language is overstated, as noted below, but the required limitation is stated. |
| Crack initiation and missing physics | 2 / 2 | Correctly states that initiation needs an additional fracture, toughness, damage, or other failure criterion absent from the data. |

**Total: 100 / 100.** Category totals: model 30/30; numerical correctness 40/40; completion and conventions 10/10; verification 15/15; scope and missing physics 5/5.

Full credit follows the specified rubric: the correct independent compatibility check satisfies the verification requirement even though an additional claimed check fails. There is no separate rubric deduction for that failed extra check. The reference script reran successfully and reproduced its recorded results; no reference correction or alternative grading is warranted.

**(a) Scientific errors**

- The direct J implementation in [submission/solve.py](submission/solve.py), lines 109–131, uses Srr, Stt, and Srt left over from the first loop's final angle, rather than recomputing stress at each integration point. It also treats polar stress components as Cartesian when forming traction. For a radial normal, the correct Cartesian traction is t_x = σ_rr cosθ − σ_rθ sinθ and t_y = σ_rr sinθ + σ_rθ cosθ. These errors invalidate that contour check, without affecting the correctly computed requested outputs or the separate compatibility check.
- Calling the yield correction “negligible” throughout the screened interval overstates the high-stress screen. At its upper endpoint σ_y/σ_e = 0.1, so the omitted plastic-factor correction is 10%. The answer's explicit statement that these screens are not finite-r error bounds is correct.

**(b) Unsupported claims**

- ANSWER.md claims direct contour values J = 0.800000 N/mm at r = 0.01, 0.05, and 0.15 mm, with agreement to 10^−11. Rerunning the unchanged supplied script produces **0.102752 N/mm at all three radii**, contradicting that claim. Its alleged confirmation of the amplitude by this particular numerical check is unsupported. See [scratch/solve.stdout](scratch/solve.stdout) and [scratch/audit.json](scratch/audit.json).

**(c) Presentation problems**

- The C_y cross-check is described as a “finite-difference derivative,” but the code directly evaluates u_y/√r; it takes no finite difference. Its numerical value is nevertheless correct.
- Calling the screened region an “annulus” is imprecise because the requested screens were evaluated only on θ = 0. The reported interval itself is correctly identified as lying on that ray.
- The discussion of an abandoned effective modulus and an alleged ambiguity in “leading field” is unnecessary; the final calculation uses the requested leading field correctly.

**Tolerance verdict:** All seven requested numbers pass: T is within 1% relative tolerance, and the other six numbers are each within 2% relative tolerance.
