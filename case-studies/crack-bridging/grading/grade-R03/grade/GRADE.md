# Grade

| rubric item | points awarded / max | reason |
|---|---|---|
| Finite geometry and coupling | 12 / 12 | Uses the original half-slit a, two connector zones, c = a + L, and the coupled piecewise pressure Σ inside the slit and Σ − t in the zones. Its finite-crack stress-intensity expression is correct; it makes no small-zone or traction-free toughness approximation. The erroneous opening relation is identified below. |
| Intrinsic front work | 7 / 7 | Correctly imposes K_tip = √(E′Γ₀), with Γ₀ = 22 J/m² at the advancing outer fronts, separately from the connector traction. Neither discards nor doubles Γ₀. |
| Elastic and local constitutive inputs | 6 / 6 | Plane-strain modulus: 3/3, E′ = E/(1 − ν²) = 141.0590278 GPa. Local law: 3/3, constant 160 MPa closing traction and a 9 µm full-gap cutoff. The kernel error is not an error in these material inputs. |
| Remote stress | 0 / 25 | Reported 158.525 MPa versus 133.730286053811 MPa: +18.5408367%, exceeding the rubric's 15% outer partial-credit limit. The bold 158.5 MPa also fails that limit. |
| Zone length per end | 0 / 25 | Reported 34.324 mm per end versus 4.08560913824553 mm: +740.1195229%, exceeding 15%. The bold 34.3 mm also fails that limit. |
| Verification | 15 / 15 | The submitted explore.py actually executes an independent weight-function quadrature for K at c = 3.4 mm, Σ = 100 MPa: analytic 2,048,805.374744487 and numerical 2,048,805.374744999 Pa√m, a relative difference of 2.50 × 10⁻¹³. The answer interprets this comparison and also checks the L = 0 limit. This satisfies the rubric's requirement of at least one meaningful independent elasticity check. The opening checks are invalid and do not establish the claimed answer accuracy; those problems are detailed below. |
| Conventions and event | 7 / 10 | Conventions: 5/5, clearly uses full relative opening, L per end, c = a + L, and the single crack-plane energy convention. Event/branch: 2/5, correctly names δ(a⁺) = δ_c and follows the tip-growth relation from initiation, but the wrong opening formula selects a physically inconsistent late event after those connectors have already failed. Its unrestricted monotonic-branch claim is also false. |

**Total: 47 / 100.**

Numerical rows use the more precise values explicitly reported in ANSWER.md. Rerunning final_solve.py gives Σ = 158.524900100573 MPa and L = 34.324303918593 mm, confirming that rounding does not explain the discrepancy.

## (a) Scientific errors

- **Incorrect finite-crack opening.** The submitted relation δ(a) = (4/E′)[Σ√(c² − a²) − t(c − a)] is not the opening of the stated two-dimensional slit under its piecewise pressure. The correct full gap is

  \[
  \delta(a)=\frac{4}{E'}\sqrt{c^2-a^2}
  \left[\Sigma-\frac{2t}{\pi}\cos^{-1}(a/c)\right]
  +\frac{8ta}{\pi E'}\ln(c/a).
  \]

  At the submitted code's root, this expression and the exact Green-kernel quadrature both give **24.3002426952 µm**, while its incorrect formula gives 9 µm. The independently integrated tip condition still gives G = 22 J/m². Thus the stress-intensity calculation is correct, but the local failure event is not.
- **Invalid opening verification.** The double integral called a general central-crack opening formula in ANSWER.md does not represent this slit elasticity problem. In explore.py, evaluating its opening integral at x = a simply integrates √(s² − a²)/√(s² − a²) = 1, reproducing the same incorrect t(c − a) term. Superposition is valid, but the supplied partially pressurized crack opening is not. For the code's c = 3.4 mm, Σ = 100 MPa check, the submitted gap is 2.292143732 µm; the correct analytic and Green calculations give 3.768343436 µm.
- **Incorrect opening saturation.** The claimed limit 4ta/E′ = 10.89 µm is false. Along the submitted tip-growth relation with Γ₀ > 0, even its own opening expression diverges as 4K₀√c/(√πE′). With Γ₀ = 0 its expression tends to 4ta(1 − 2/π)/E′, about 3.957 µm. The correct finite-crack opening also does not saturate at the claimed value.
- **Incorrect global load monotonicity.** Its correct stress relation reaches a maximum Σ ≈ 161.010172 MPa at L ≈ 239.620248 mm, then decreases toward 160 MPa from above. It cannot increase for all L while tending to 160 MPa from above. The load does increase up to the reported root; the physical inconsistency there is that the original connectors reached their cutoff earlier, at the reference event.

## (b) Unsupported claims

- The opening calculations do not agree to machine precision as claimed. Rerunning verify_super.py emits an IntegrationWarning and reports 9.0003134704 µm against 8.9999998929 µm, a relative discrepancy of approximately 3.48 × 10⁻⁵. Both calculations use the wrong opening mechanics.
- Perturbing δ_c by ±0.3% does reproduce small changes in the submitted model's outputs, but sensitivity to an input does not demonstrate accuracy against the reference. The claim that both answers are consequently within ±0.3% is contradicted by their actual errors.
- The claimed independent validation and uniqueness of the physical failure event are not established by solving residuals that reuse the erroneous opening relation.

## (c) Presentation problems

- The requested quantities, full-gap convention, and per-end length are clearly labeled; there is no major formatting problem. The supporting K = 1.7616181e6 in the checks omits its units, Pa√m.
- The checks mix valid stress-intensity evidence with invalid opening evidence without distinguishing what each establishes.

The reference implementation was rerun in scratch and reproduced both recorded targets while passing all its assertions. The analytic opening and Green quadrature also agree at the submitted root. No mechanical error in the reference was found, so no alternate corrected-reference grade is needed. Rerun outputs and the local audit are retained in scratch/.

**Tolerance verdict:** Neither requested number is within ±0.3%; acceptable intervals are Σ = 133.3290952–134.1314769 MPa and L = 4.07335231–4.09786597 mm per end.
