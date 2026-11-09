# Grade

| rubric item | points awarded / max | reason |
|---|---:|---|
| Leading tensor compliance and transverse strain | 10 / 10 | Correctly obtains effective constants T and ν* = 1/2 − (1/2 − ν)T/E. This preserves (1 − 2ν*)/T = (1 − 2ν)/E, hence the elastic volume response, and permits nonzero ε_zz under plane stress. |
| Traction-free tensile leading field | 10 / 10 | Correct mode I field, r^−1/2 stresses, r^1/2 displacements, and polar angular coefficients after the constitutive reduction. |
| J normalization | 5 / 5 | Uses J = K²/T and K = √(JT); direct contour integration confirms this normalization. |
| Inversion of transverse measurement | 5 / 5 | Uses C_zz = −2ν*√(J/(2πT)) and solves for the unique admissible T. The signed monotonicity statement has a prose error noted below, but the inverse equation and solution are correct. |
| T | 12 / 12 | Reports 4816.47 MPa versus 4816.472187302739 MPa: relative error 0.0000454%, within 1%. |
| C_rr at 60° | 5 / 5 | Reports +0.0039420 mm^1/2 versus +0.003942048237733179: relative error 0.001224%, within 2%. |
| C_θθ at 60° | 5 / 5 | Reports +0.00063318 mm^1/2 versus +0.0006331781350710297: relative error 0.0002945%, within 2%. |
| C_γrθ at 60° | 5 / 5 | Reports +0.0057311 mm^1/2 versus +0.005731131133456488: relative error 0.0005432%, within 2%; correctly uses engineering shear 2ε_rθ. |
| C_y at 90° | 5 / 5 | Reports +0.0091390 mm^1/2 versus +0.009139018941480233: relative error 0.0002073%, within 2%. |
| Lower radius endpoint | 4 / 4 | Reports 0.01000 mm, agreeing to rounding. All three leading principal-strain coefficients are (0.00264150818, 0.00264150818, −0.0050000) mm^1/2, so ε_zz controls the lower bound. Explicitly identifies a nonempty interval. |
| Upper radius endpoint | 4 / 4 | Reports 0.06133 mm versus 0.06132522854990914 mm: relative error 0.007781%, within 2%. Correct stress-invariant screen; the 0.200 mm cap is inactive. |
| Completion | 5 / 5 | All seven requested numerical outputs are identifiable, including both interval endpoints. |
| Units and conventions | 5 / 5 | Correct MPa, mm, and mm^1/2 units for the requested results (2); tensor versus engineering shear convention (2); polar and displacement signs (1). |
| Verification: independent condition | 5 / 5 | Computes J directly from W, tractions, and Cartesian finite differences of the displacement field, rather than substituting into the measured C_zz. |
| Verification: quantitative evaluation | 5 / 5 | Rerun of the submitted contour calculation gives J = 0.7999999999966, 0.8000000000111, 0.8000000000107, and 0.8000000000111 N/mm at R = 0.02, 0.05, 0.10, and 0.20 mm. |
| Verification: interpretation | 5 / 5 | Correctly interprets contour agreement with the supplied J as validating the amplitude, angular/displacement field, and J–K relation together. |
| Scope of screening approximation | 3 / 3 | Explicitly states that these are defined checks of the leading field, not error estimates for the full finite-r solution, and discusses kinematic and matching limits. Some accompanying accuracy statements are overstated, as noted below. |
| Missing initiation physics | 2 / 2 | Correctly explains that initiation needs an additional fracture resistance or failure/damage criterion. |

**Total: 100 / 100.**

The reference rerun passed all its assertions, and no reference error was found. The issues below do not change the requested leading-field results or invalidate a scored rubric criterion; no additional deductions outside the rubric are imposed. Local numerical evidence is saved in `scratch/audit_results.json`.

**(a) Scientific errors**

- The signed calibration function C_zz(T) is strictly increasing on 0 < T < E; its magnitude decreases. The answer incorrectly calls the signed left-hand side strictly decreasing. Uniqueness and the computed root remain correct.
- The auxiliary remote-loading assertions J = πσ∞² and K_eff = σ∞√(πT) are dimensionally inconsistent with the specified units and lack a geometry factor. They are not used in the correct local normalization J = K²/T.
- The claim of about 2% finite-r correction at r = 0.01 mm is contradicted by the submitted constitutive-substitution calculation: √r ε_zz = −0.004806664875, a 3.8667% difference from −0.005000 relative to the leading coefficient. At r = 0.20 mm the analogous difference is 17.2924%. These are substitutions of leading stresses into the complete law, not a complete finite-r field solution. Negligible matching corrections do not eliminate the constitutive intercept correction.
- Absolute strain corrections from the dropped σ_y/σ_e term are O(1), while their size relative to the singular leading strain is O(√r). The answer's wording conflates these orders, although it correctly identifies bounded corrections and unchanged leading coefficients.

**(b) Unsupported claims**

- The stated agreement of two executed T solvers is not demonstrated by the supplied code: it executes Brent's method, but not a separate quadratic inversion. The written quadratic is correct, and this extra claim is unnecessary for verification credit.
- Claims of a finite-r accuracy window or “self-consistent description” beyond the defined screens are not established by the calculations. Substituting leading stresses into the full constitutive law does not establish equilibrium and compatibility of a full finite-r solution. The answer does also explicitly acknowledge that the screens are not full-solution error estimates.

**(c) Presentation problems**

- The code prints full-law strain diagnostics at the screening endpoints without clearly distinguishing them from the leading-field strains used to determine those endpoints. For example, it prints 0.04807 at the lower endpoint, whereas the requested leading strain screen is exactly 0.05000 there. The reported endpoint itself is correct.

**Tolerance verdict:** All seven requested numbers pass: T is within 1% and the other six are within 2%; the reported interval is nonempty.
