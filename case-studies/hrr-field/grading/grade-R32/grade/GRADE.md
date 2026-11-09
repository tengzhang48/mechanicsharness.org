# Grade

| rubric item | points awarded / max | reason |
| --- | ---: | --- |
| Constitutive model and inner limit: hardening dominance | 10 / 20 | Correctly selects the nonlinear inner field, but only asserts that it dominates (ANSWER.md, lines 21–24). This earns the rubric's explicit “correct choice with little reasoning” score. The missing explanation is that, for positive p and n = 4, elastic strain scales as r^−p while hardening strain scales as r^−4p, so their ratio grows as r^−3p. |
| Constitutive model and inner limit: plastic plane strain and effective stress | 10 / 10 | `notch.py` implements Q² = 3[f″ + s(2−s)f]²/4 + 3(s−1)²(f′)² and the correct plane-strain compatibility equation. This is the effective stress associated with σzz = (σrr + σθθ)/2 in the leading plastic limit; the constitutive factors are correct. |
| Geometry: free faces | 10 / 10 | Correct material interval and free-face conditions f(±2π/3) = f′(±2π/3) = 0. |
| Geometry: symmetry and amplitude normalization | 5 / 5 | Uses f′(0) = f‴(0) = 0, f(0) = 1, and a positive opening amplitude, with consistent Airy stresses. |
| Geometry: leading tensile branch and finite energy | 5 / 5 | `eig.py` selects the opening branch by continuation from the crack. The answer checks p < 2/(n+1) = 0.4. The reproduced solution has positive f throughout the sampled interior and radial energy power 2−5p = 1.219123698570 > 0. |
| Numerical exponent p | 15 / 15 | Rounded p = 0.15618 has absolute error 0.000004739714 < 0.0003 against 0.156175260286. Its parenthesized value matches every printed reference digit, and the code reproduces it. |
| Numerical stress ratio R | 15 / 15 | Rounded R = 0.3733 has absolute error 0.000037686688 < 0.002 against 0.373262313312. The formula R = [s + f″(0)]/[s(s−1)] is correct, and the code reproduces the parenthesized value. |
| Quantitative verification: independent method or analytic/limiting check | 10 / 10 | The executed crack-limit checks reproduce s = (2n+1)/(n+1) for n = 3, 4, 5, 9, with maximum error 5.82×10^−14. This successful analytic limiting check earns full credit. RK4 versus adaptive integration also reproduces Δs = 1.18×10^−14 and Δf″(0) = 7.55×10^−15. |
| Quantitative verification: residuals or refinement | 5 / 5 | `bccheck.py` reproduces f(β) = 7.022×10^−13 and f′(β) = 4.278×10^−13. Quantified refinement is also successful, although its claimed precision is overstated below. |
| Justification and scope | 3 / 5 | Identifies the monotonic proportional total-deformation setting and correctly explains amplitude cancellation. It does not explain that global matching to the remote loading and specimen determines the unrequested amplitude K. Award the rubric's one-component score. |

**Total: 88 / 100.** The deductions are 10 points for the incomplete dominance explanation and 2 points for the missing global-matching explanation. The minor verification and presentation issues below incur no additional deductions under this rubric.

**(a) Scientific errors**

- Calling the nonlinear strain term “creep” is incorrect terminology for the supplied rate-independent total-deformation law. No time or strain-rate constitutive equation is specified. The submitted numerical model nevertheless uses the correct hardening field.
- No error was found in the implemented plane-strain equation, stress convention, face conditions, or requested numerical results.

**(b) Unsupported claims**

- The claimed change in s of less than 10^−12 between 1,000 and 8,000 RK4 steps is contradicted by rerunning those solves: the change is 3.028688411177427×10^−12. This remains excellent convergence for the requested tolerance.
- Crack-limit agreement “to 1e-14” for all four listed n values is overstated: rerun errors range from 2.02×10^−14 to 5.82×10^−14. The analytic check still succeeds comfortably.
- The vectorized alternate-route description says `brentq` solves the inner condition f(β) = 0. In `vnotch1d_fast.py`, that inner root is linearly interpolated; `brentq` solves the outer condition. The separate scalar `vnotch1d.py` does implement an inner Brent solve. The fast route's reported s = 1.8438247424 is reproduced and agrees with the principal result to about 2.7×10^−9.

**(c) Presentation problems**

- The answer reuses α for the angular boundary 2π/3, although the problem already uses α = 0.02 for the material coefficient. Its material interval is correct; a separate angular symbol would avoid ambiguity.
- `bccheck.py` prints the shear-stress formula with f(α) where it should say f′(α). The actual calculation correctly uses f′(α).

**Reference and reproduction:** The copied `reference.py` runs successfully and reproduces its recorded results. No reference error was found, so no alternative corrected-reference grade is needed. Submission reruns, the reference rerun, and a precise refinement/branch audit are retained in [scratch/](scratch/), including [audit_results.json](scratch/audit_results.json) and the individual run logs. Neither the original reference nor the original submission was modified.

**Tolerance verdict:** Both requested numbers are within their stated absolute tolerances, including the rounded headline values p = 0.15618 and R = 0.3733.
