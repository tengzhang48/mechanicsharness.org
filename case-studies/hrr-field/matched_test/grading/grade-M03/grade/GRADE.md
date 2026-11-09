# Grade

| rubric item | points awarded / max | reason |
|---|---:|---|
| Leading tensor compliance | 10 / 10 | Correctly preserves the elastic bulk response, modifies the deviatoric compliance, and obtains E_eff = T and ν_eff = 0.486238650893, equivalent to ν_t = 1/2 − (1/2 − ν)T/E. Allows nonzero ε_zz. |
| Tensile traction-free leading field | 10 / 10 | Correct plane-stress Williams mode-I field, angular stresses, r^−1/2 stresses, and r^1/2 displacements, justified after reducing the constitutive law. |
| J normalization | 5 / 5 | Correct J = K²/T; the submitted direct energy-flux calculation independently supports this normalization. |
| Transverse inverse measurement | 5 / 5 | Correct ε_zz = −A_tr(σ_xx + σ_yy), inverse equation, and unique admissible root with 0 < T < E; the submitted bracketing solve reproduces it. |
| T | 12 / 12 | Reports 4816.47 MPa versus 4816.472187302739 MPa; relative error 0.0000454%, within 1%. |
| C_rr(60°) | 5 / 5 | Reports +0.0039420 mm^1/2 versus +0.003942048237733179; relative error 0.001224%, within 2%. |
| C_θθ(60°) | 5 / 5 | Reports +0.00063318 mm^1/2 versus +0.000633178135071030; relative error 0.0002945%, within 2%. |
| C_γrθ(60°) | 5 / 5 | Reports +0.0057311 mm^1/2 versus +0.005731131133456488; relative error 0.0005432%, within 2%. Correct engineering shear factor of two. |
| C_y(90°) | 5 / 5 | Reports +0.0091390 mm^1/2 versus +0.009139018941480233; relative error 0.0002073%, within 2%. |
| Lower radius endpoint | 4 / 4 | Reports 0.01000 mm, matching the reference. Examines all three principal strains, identifies ε_zz as binding, and explicitly reports a nonempty interval. The code correctly squares the coefficient-to-strain ratio despite a written typo. |
| Upper radius endpoint | 4 / 4 | Reports 0.061325 mm versus 0.06132522854990914 mm; relative error 0.0003727%, within 2%. Correct effective-stress bound; the 0.200 mm cap is nonbinding. |
| Completion | 5 / 5 | All seven requested outputs, including both interval endpoints, are identifiable. |
| Units and conventions | 5 / 5 | Correct units (2), tensor versus engineering shear (2), and polar/displacement signs (1). |
| Verification: independent condition | 5 / 5 | verify_J.py evaluates W and displacement derivatives in the contour integral instead of substituting into C_zz. solve.py also reconstructs displacement from the leading strain field. |
| Verification: quantitative evaluation | 5 / 5 | Unmodified verify_J.py gives plane-stress J E/K² = 0.999999660–0.999999766 at four Poisson ratios. Its submitted flux function, rerun at the calibrated parameters and radii 0.001–0.050 mm, gives J ≈ 0.799999721 N/mm. solve.py gives C_y = 0.009139018941480233 and reconstructed C_y = 0.009139018941480228. |
| Verification: interpretation | 5 / 5 | Correctly interprets agreement as validating J normalization and the displacement coefficient, independently of the transverse calibration. One successful check suffices under the rubric. |
| Approximation scope | 3 / 3 | Explicitly treats the interval as defined checks on the leading field, not error bars for the full finite-r solution. The additional accuracy wording is too strong, as noted below, but the required distinction is present. |
| Missing initiation physics | 2 / 2 | Correctly states that initiation requires a fracture, damage, or other failure criterion; constitutive parameters and J alone do not supply a threshold. |

**Total: 100 / 100.** Category totals: model 30/30, numerical correctness 40/40, completion/conventions 10/10, verification 15/15, scope/missing physics 5/5.

The reference code reran successfully and reproduced the recorded results and checks. No reference correction is warranted. The observations below do not remove points under the specified rubric: the requested leading-field outputs and the relevant computations are correct.

Scientific errors:

- The written lower-bound step says r ≥ 0.005/0.05 = 0.010000 mm. It should say r ≥ (0.005/0.05)² = 0.010000 mm. The unsquared ratio is 0.1 mm^1/2, not a radius. The submitted code uses the correct square.

Unsupported claims:

- The statements that the reduction is “accurate” or “good to ~10%” overreach the stress screen. That screen bounds the omitted bracket correction, not the relative error of every total-strain component or the full finite-r solution. For example, at σ_xx = σ_yy = 100 MPa, σ_zz = 0, the leading law gives ε_zz = −0.0201906554, while the complete law at the same stress gives −0.0182573041, a 10.59% difference relative to the latter. This comparison is not a finite-r boundary-value solution. The answer does separately make the required statement that the screens are not full-solution error bars.
- The claimed measurement-dominated uncertainty of approximately 0.1% in T is unsupported: the problem gives no uncertainty model for C_zz. Its displayed value, −0.005000, also has four significant figures, not the claimed five.

Presentation problems:

- K denotes the stress-intensity amplitude and is also used for the bulk modulus in 1/E_eff = 1/(3G_eff) + 1/(9K). Distinct symbols would avoid ambiguity.
- The auxiliary in-plane strain coefficient is printed as 0.0026408; the submitted code gives 0.00264150817794 mm^1/2. This minor discrepancy does not alter the binding out-of-plane strain screen or any requested output.

Tolerance verdict: all seven requested numbers meet their stated relative tolerances (1% for T, 2% for the other six); the reported closed interval [0.01000, 0.061325] mm is nonempty.
