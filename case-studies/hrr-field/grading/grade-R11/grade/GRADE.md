# Grade

| rubric item | points awarded / max | reason |
| --- | --- | --- |
| Stress exponent: numerical value | 10 / 10 | Reported p = 0.133333 (also 0.1333333333); even the coarser value differs from the key by 3.33333×10⁻⁷, below 10⁻⁵. |
| Stress exponent: contour-flux scaling | 2 / 5 | Correctly states the stress/displacement powers and W ∝ r⁻¹. However, p is justified through the angular eigenvalue, without explaining why nonzero, radius-independent J forces 1 − (n+1)p = 0. Award the rubric's incomplete-scaling credit. |
| Mechanics: equilibrium and compatibility | 10 / 10 | Correct Airy stress components enforce equilibrium; describes compatibility shooting and supplies a consistent strain–displacement reconstruction. Reported solved angular quantities I = 4.818379 and normalized hoop amplitude 2.335257 agree with the key. |
| Mechanics: plane-strain closure and constitutive factors | 10 / 10 | Correct σzz = (σrr+σθθ)/2 and σe² = 3(σrr−σθθ)²/4 + 3σrθ². Reconstruction uses normal-strain factor 3/4 and tensor-shear factor 3/2, with no added elasticity. |
| Mechanics: symmetry, free faces, normalization | 5 / 5 | States f′(0) = f‴(0) = 0 and f(π) = f′(π) = 0, and explicitly normalizes max σ̃e = 1. |
| Tensile coefficient H: numerical value | 20 / 20 | Reported H = 1.8936 differs from 1.893569725213 by 0.000030274787, below 0.005. |
| Tensile coefficient H: J calibration | 5 / 5 | Uses the full-circle contour and H = σ̃θθ(0)/I^(1/(n+1)); the dimensional reconstruction consistently uses αε₀σ₀ and the given J convention. |
| Energy coefficient C: numerical value | 15 / 15 | Reported C = 0.4573 differs from 0.457305783974 by 0.000005783974, below 0.002. |
| Energy coefficient C: radial integration and thickness | 5 / 5 | Correctly integrates W ∝ r⁻¹ with area element r dr dθ to obtain U(R) proportional to R, and uses energy per unit thickness in C = U/(JR). |
| Verification: meaningful second calculation | 10 / 10 | Reports full-circle dimensional J from a separately reconstructed displacement and both analytic and finite-difference Cartesian gradients, with J/Jgiven near unity. The submitted reconstruction formulas support this check; the diagnostic described below corroborates them. |
| Verification: equation residuals or convergence | 5 / 5 | Reports a 1.5×10⁻⁹ discrepancy between the two angular-displacement constructions, with scale about 2.2. The final verification script explicitly calculates this shear-compatibility residual; the diagnostic reproduces 1.517×10⁻⁹. |

**Total: 97 / 100.**

The reference was rerun in `scratch/` and reproduced its recorded results and checks. No reference error was found, so no alternate corrected-reference grade is needed.

Both unmodified submission reruns (`compute.py` and `verify3.py`, copied into `scratch/`) fail at `import hrr`. To audit the supplied reconstruction formulas, `scratch/audit_submission_checks.py` supplies angular fields from the reference and omits the submission's eigenvalue shooting. This diagnostic obtains H = 1.89356973, C = 0.45730578, analytic J/Jgiven = 0.9999999999, and maximum analytic-versus-FD gradient difference 4.294×10⁻¹⁷. It corroborates the downstream mechanics, but does **not** establish that the missing angular solver or the original scripts were successfully rerun. Missing supplementary material alone receives no deduction; verification credit rests on the concrete numerical checks reported in the answer and the supplied reconstruction formulas.

Scientific errors:

- The stated equality s = (2n+1)/(n+1) = **13/7.5** is arithmetically wrong: it should be **14/7.5**. The stated decimal s = 1.86667, formula, and final p are correct, so this typo does not change the numerical credit.
- No constitutive, plane-strain invariant, shear-factor, J-convention, or final energy-integration error was found.

Unsupported claims:

- The claim that the angular boundary conditions admit a nontrivial solution **only** at the stated s is not established by searching a finite interval and selecting the first shooting bracket. The relevant leading-field exponent is nevertheless correct.
- The original solver tolerance, the eigenvalue difference of 5×10⁻¹⁵, and the asserted 2001→32001 grid-convergence history cannot be independently reproduced from the supplied files. There is no submitted convergence loop or sequence of outputs. These claims are not needed for the verification points awarded above.

Presentation problems:

- The submitted scripts lack their required `hrr` dependency and therefore are not independently runnable from this folder.
- The angular compatibility ODE and the integrand defining I are abbreviated in the prose, reducing the method's self-containedness.

Requested-number tolerance verdict: **all three pass** — |Δp| = 3.33333×10⁻⁷ < 10⁻⁵; |ΔH| = 3.02748×10⁻⁵ < 0.005; |ΔC| = 5.783974×10⁻⁶ < 0.002.
