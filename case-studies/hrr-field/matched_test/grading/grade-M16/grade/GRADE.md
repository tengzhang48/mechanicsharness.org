# Grade

| rubric item | points awarded / max | reason |
|---|---:|---|
| Leading tensor compliance and transverse strain | 10 / 10 | Correct effective pair E* = T and nu* = 1/2 - (1/2 - nu)T/E; retains elastic volume response and allows nonzero epsilon_zz. The transverse relation includes both elastic and deviatoric contributions. |
| Traction-free tensile leading field | 10 / 10 | Correct plane-stress mode-I Williams field after the constitutive reduction, with stresses proportional to r^(-1/2) and displacements proportional to sqrt(r); explains the finite-J exponent. Angular stresses and rotations in the code are correct. |
| J normalization | 5 / 5 | Uses J = K_loc^2/T, giving K_loc = 62.073970 MPa sqrt(mm). The submitted direct flux calculation independently supports this normalization. |
| Inversion using C_zz | 5 / 5 | Correct equation C_zz(0) = -2 nu* sqrt(J/(2 pi T)) and successful admissible root T = 4816.472187 MPa. Its magnitude decreases strictly for 0 < T < E, so the calibrated root is unique. |
| Numerical T | 12 / 12 | Reports 4816.5 MPa (also 4816.47); relative error using 4816.5 is 0.000577%, below 1%. |
| Numerical C_rr(60 degrees) | 5 / 5 | Reports +0.003942 mm^(1/2); relative error 0.001224%, below 2%. |
| Numerical C_thetatheta(60 degrees) | 5 / 5 | Reports +0.000633 mm^(1/2); relative error 0.02813%, below 2%. |
| Numerical C_gamma_rtheta(60 degrees) | 5 / 5 | Reports +0.005731 mm^(1/2) for engineering shear; relative error 0.002288%, below 2%. |
| Numerical C_y(90 degrees) | 5 / 5 | Reports +0.009139 mm^(1/2); relative error 0.0002073%, below 2%. Correctly uses u_y = u_r on this ray and C_y = 2 C_rr there. |
| Numerical lower radius endpoint | 4 / 4 | Reports 0.0100 mm, matching the reference to roundoff. Explicitly examines all three principal strains and identifies epsilon_zz as controlling. Interval is explicitly nonempty. |
| Numerical upper radius endpoint | 4 / 4 | Reports 0.0613 mm; relative error 0.04114%, below 2%. Correct stress-invariant screen and minimum with 0.200 mm; explicitly reports the closed nonempty interval. |
| All requested outputs identifiable | 5 / 5 | All seven numerical outputs are reported, including both interval endpoints. |
| Units | 2 / 2 | Correct MPa, mm, and mm^(1/2) units. |
| Tensor versus engineering shear | 2 / 2 | Defines C_gamma_rtheta using 2 epsilon_rtheta and applies the corresponding factor of two. |
| Polar and displacement signs | 1 / 1 | Consistent basis rotations, positive requested coefficients, tensile amplitude, and displacement convention. |
| Verification: independent condition | 5 / 5 | Directly evaluates the J contour integral from energy density, traction, and finite-difference displacement derivatives; this is independent of substitution into the C_zz calibration. |
| Verification: quantitative evaluation | 5 / 5 | Rerunning the submitted physical-unit calculation gives J = 0.7999999934947 N/mm at radii 0.005, 0.020, 0.050, and 0.100 mm (maximum relative error 8.13e-9). The complete-law uniaxial check also reproduces its bilinear target. |
| Verification: interpretation | 5 / 5 | Explicitly compares the contour values with the supplied 0.800 N/mm and uses the agreement to verify normalization. Also explains the separate displacement-coefficient consistency check. |
| Scope of the screens | 3 / 3 | Recognizes the specified leading-field working interval and explicitly says the screens are not finite-r solution error estimates or a guarantee of continuum validity at the tip. |
| Missing initiation physics | 2 / 2 | Correctly states that initiation needs an additional fracture resistance, critical separation, or damage criterion; field amplitude and the strain screen alone do not supply it. |

**Total: 100 / 100.** The reference rerun reproduces its recorded results and passes its independent assertions; no reference correction or alternative score is warranted.

Scientific errors:

- No substantive errors in the requested leading-field mechanics, numerical results, or screening calculation.
- The explanatory sentence claiming to derive E* = T from the elastic trace and the identity (1 + nu*)/E* - nu*/E* = 1/E* is incomplete. The trace alone fixes (1 - 2 nu*)/E*; the original deviatoric compliance is also needed to determine E*. The correct effective pair is explicitly given and used, satisfying the rubric's full-credit reduction criterion.

Unsupported claims:

- None needed for the awarded mechanics credit. The numerical outputs, full-law comparison at prescribed leading stresses, uniaxial reduction, and independent physical-unit J check are supported by the submitted calculations and reruns.

Presentation problems:

- Section 3(i)'s asymptotic wording is imprecise: the omitted yield-intercept contribution is O(1) in strain and O(sqrt(r)) in sqrt(r)-normalized strain. The earlier leading-limit explanation is correct.
- Calling the screened range an "annulus" is loose wording: the requested and demonstrated screening interval is on theta = 0. The answer's explicit interval and calculations correctly identify that ray.

Tolerance verdict: **All seven requested numbers pass: T is within 1%, and all six remaining outputs are within 2%.**
