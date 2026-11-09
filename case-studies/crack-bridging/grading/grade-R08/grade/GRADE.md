| rubric item | points awarded / max | reason |
|---|---:|---|
| Finite geometry and coupling | 12 / 12 | Uses the finite symmetric crack with half-length c = a + L and constant closing traction on a < |x| < c. The coupled opening and tip-intensity formulas match the reference; no small-zone approximation is used. |
| Intrinsic front work | 7 / 7 | Enforces K_tip = sqrt(E′ Γ₀), with Γ₀ = 22 J/m² at each advancing tip, separately from connector work. No doubling or replacement by total plateau toughness occurs. |
| Elastic and local constitutive inputs | 6 / 6 | Modulus: 3 / 3 for E′ = E/(1 − ν²) = 141.059 GPa. Traction/opening: 3 / 3 for constant 160 MPa closure and the full relative opening, with cutoff 9 μm and connector work 1440 J/m². |
| Remote stress | 25 / 25 | Reported 133.73 MPa versus 133.730286053811 MPa: absolute relative error 0.0002139035%, within ±0.3%. The copied closed-form solver reproduces 133.730286 MPa. |
| Zone length per end | 25 / 25 | Reported 4.086 mm per end versus 4.08560913824553 mm: absolute relative error 0.0095667926%, within ±0.3%. The copied solver reproduces 4.085609 mm. |
| Verification | 15 / 15 | Executes and interprets the known Γ₀ = 0 limit, reproducing S/σ₀ = 0.823969557 with a 1.11 × 10⁻¹⁶ difference between formulations. A separate finite-crack spectral elasticity calculation also reproduces the quoted grid-dependent bracket. This exceeds merely checking imposed tip/opening residuals. |
| Conventions and event | 10 / 10 | Conventions: 5 / 5 for full gap, per-end L = (η − 1)a, and crack-plane work counted once. Event: 5 / 5 for imposing δ(a) = δ_c at the trailing connector and K_tip = K₀ at c, on the physically consistent monotonic growth branch. The stated formulas confirm no earlier load maximum. |

**Total: 100 / 100.**

Reran copies of `closed.py`, `verify.py`, `mono.py`, and `final.py` inside `scratch/submission/`. The spectral runs give (N, Σ in MPa, L in mm): (300, 133.82971, 4.132229), (600, 133.60406, 4.051940), and (1000, 133.89042, 4.117758), matching the submission's recorded results. Those coarse lengths are a check, not the reported final answers. The reference rerun passes its assertions and independently recovers the targets by Green-function quadrature; no reference correction or alternate grade is warranted.

**(a) Scientific errors**

- No scientific error affects the final equations, implementation of the reported solution, or requested numbers.
- The parenthetical claim that the “K-concept” does not apply is too broad. The original-flaw small-zone approximation is inappropriate, but the local tip intensity and K_tip²/E′ = Γ₀ remain applicable, as the submission itself correctly uses.
- The initial docstring in `solve.py` incorrectly places first connector failure at |x| = c with δ(c) = δ_c. The opening vanishes at that tip. The final answer and the routines actually used for the result correctly evaluate the cutoff at |x| = a, so this obsolete comment causes no rubric deduction.

**(b) Unsupported claims**

- The spectral scatter is attributed specifically to the logarithmic junction singularity without an error analysis isolating that cause. The rerun supports the numerical bracket, but the three grids do not demonstrate convergence to all reported digits.
- `mono.py` checks discrete states beginning at η = 1.2; those samples alone do not prove monotonicity continuously from initiation at η = 1. The displayed analytic load relation does corroborate the claim: its derivative stays positive through the cutoff, where d(Σ/σ₀)/d ln η = 0.2150159321. No event credit depends on an incorrect branch.

**(c) Presentation problems**

- Exploratory scripts retain inconsistent root/tip labels, particularly `solve.py`; the final answer resolves the event and units clearly.
- The two nonlinear solves use the same elasticity equations, so their agreement checks numerical consistency rather than providing independent elasticity validation. The separately shown spectral calculation and zero-front-work limit supply the substantive checks.

**Tolerance:** Both requested final numbers are within ±0.3% of the unrounded reference targets.
