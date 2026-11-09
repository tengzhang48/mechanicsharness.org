# Grade

| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: high-stress tensor compliance | 10 / 10 | Correct spherical and deviatoric reduction, with E_eff = T and ν_eff = 1/2 − (1/2 − ν)T/E. The elastic bulk response is retained, and ε_zz is correctly nonzero. The stated A_hyd = 1/(9K_bulk) is correct because S is the stress trace. |
| Model: traction-free tensile leading field | 10 / 10 | Correct plane-stress mode-I field, r^−1/2 stresses and strains, and r^1/2 displacements. The harmonic polar stress functions equal the reference functions and satisfy traction-free faces. |
| Model: normalization using J | 5 / 5 | Correctly distinguishes K_A from the standard stress-intensity factor and uses J = 2πK_A²/T = K_loc²/T. |
| Model: transverse inverse calibration | 5 / 5 | Uses C_zz(0) = −2ν_eff√[J/(2πT)] and solves for T = 4816.472 MPa. The equation has a unique positive root, and this root lies in 0 < T < E. |
| Numerical: T | 12 / 12 | Reports 4816.47 MPa versus 4816.472187302739 MPa; relative error 0.0000454%, within 1%. |
| Numerical: C_rr at 60° | 5 / 5 | Reports +0.0039420 mm^1/2 versus +0.003942048237733179; relative error 0.001224%, within 2%. |
| Numerical: C_θθ at 60° | 5 / 5 | Reports +0.00063318 mm^1/2 versus +0.0006331781350710297; relative error 0.000295%, within 2%. |
| Numerical: C_γrθ at 60° | 5 / 5 | Reports +0.0057311 mm^1/2 versus +0.005731131133456488; relative error 0.000543%, within 2%. Correct engineering shear and sign. |
| Numerical: C_y at 90° | 5 / 5 | Reports +0.0091390 mm^1/2 versus +0.009139018941480233; relative error 0.000207%, within 2%. The stated radial displacement formula gives the correct u_y on this ray. |
| Numerical: lower radius endpoint | 4 / 4 | Reports 0.010000 mm, agreeing with the reference. Examines all three principal strains and correctly identifies the out-of-plane strain as controlling. |
| Numerical: upper radius endpoint | 4 / 4 | Reports 0.061325 mm versus 0.06132522854990914 mm; relative error 0.000373%, within 2%. Correct stress-invariant screen, geometric cap, and explicitly nonempty interval. |
| Completion: all seven outputs | 5 / 5 | T, the three angular strain coefficients, the displacement coefficient, and both radius endpoints are all reported and identifiable. |
| Conventions: units | 2 / 2 | Correct MPa, mm^1/2, and mm units for the requested outputs. |
| Conventions: tensor and engineering shear | 2 / 2 | Explicitly defines engineering shear as twice tensor shear and applies the factor correctly. |
| Conventions: polar and displacement signs | 1 / 1 | Positive requested coefficients and the polar basis/displacement signs agree with the reference. |
| Verification: independent condition | 5 / 5 | Section 7.1 specifies direct evaluation of the defined J contour integral, an independent flux check rather than reuse of the measured C_zz equation. |
| Verification: quantitative evaluation | 5 / 5 | Reports J = 0.8000000000 N/mm at four radii. Grader recalculation gives 0.800 N/mm within 3 × 10^−16. The separate full-law finite-radius C_zz sequence in §7.3 is also reproduced. |
| Verification: interpretation | 5 / 5 | Explains that the direct integral confirms amplitude normalization and that the full-law sequence confirms the yield contribution is subleading. The unauditable additional check residuals are not needed for this credit. |
| Scope: meaning of the screens | 3 / 3 | Correctly describes the interval as simultaneous checks of the stated leading field and expressly says these are not finite-radius solution error estimates. |
| Missing physics: crack initiation | 2 / 2 | Correctly states that initiation requires an additional fracture/failure criterion, such as critical J, toughness, or a damage law. |

**Total: 100 / 100.** Category totals: model 30/30; numerical correctness 40/40; completion and conventions 10/10; verification 15/15; scope and missing physics 5/5.

Scientific errors:

- No substantive scientific error affecting the requested solution was found. The leading compliance, field, normalization, inversion, displacement, and three-principal-strain screen are correct.

Unsupported claims:

- The detailed execution and residual claims in §7.2, including the reported equilibrium and compatibility checks, cannot be audited from the submitted files: only `ANSWER.md` is present. No credit depends on that block. The independently specified and quantitatively reported J check, corroborated by recalculation, supplies the required verification; submitting code is not itself a rubric requirement.

Presentation problems:

- In §4, the expressions written as ε_rr, ε_θθ, and ε_rθ using Σ should instead be labeled strain coefficients, or should include r^−1/2. The following C formulas and requested numerical values are correct.
- The auxiliary κ_eff is printed as 1.69137; the stated formula gives 1.69135781, which rounds to 1.69136. This minor discrepancy does not affect any requested value or convention.

Reference assessment: the reference script passed its assertions, and an independent calculation using the submission's harmonic stress functions and the original leading tensor compliance reproduced all seven reference values. No reference correction or alternate score is warranted. Local audit details are in [scratch/audit_results.json](scratch/audit_results.json).

Tolerance verdict: **all seven requested numbers are within tolerance** (T within 1%; the other six within 2%), with correct signs and an explicitly nonempty screening interval.
