# Grading profile

All five requested numerical results pass the key's bands and reproduce exactly. Part (a) is partly right because its additional description gives an incorrect initial rate and an incorrectly normalized late-time traction. The final numerical method is appropriate; several explanations and archived experiments contain errors.

## 1. Answers

| Question | Judgment | Reason |
|---|---|---|
| (a) | **Partly right** | The modal decay, first-mode shape, λ₁ = 11.7426678191 and λ₂ = 264.634203163 are correct within the key's exact/2% bands; the asserted zero initial interior rate and the supplied late-time stress amplitude are wrong (D1–D2 below). |
| (b) | **Right** | The limiting boundary is a traction-free, crack-like diffusion wedge; K∞/(σ₀√(πh)) = 1.12152225660 is inside 1.09–1.15. |
| (c) | **Right** | The slowest-mode e-folding time is 1/λ₁ = 0.0851595238326 t₀, inside 0.05–0.15 and exactly consistent with the solver's λ₁. |

## 2. Supporting state

**Right for the state requested by the problem.** The normalized first-mode midpoint value, 0.8019565042704974, agrees with the key's 0.8020 within 2%. Its reported quarter/three-quarter values and monotone shape also agree with the delivered eigenvector. The dimensionless evolution equation σ_s = K[σ_xxx], time normalization, and three traction boundary conditions are correct. The incorrect extra higher-mode/amplitude claims are catalogued separately.

## 3. Method

**Right: singular-integral Chebyshev collocation with resolvent/SVD pole finding.** The final code represents f''' with the appropriate inverse-square-root tip weight, handles the PV term analytically with polynomial subtraction, evaluates the regular image terms numerically, and builds f(0) = f'(1) = f''(1) = 0 into the triple antiderivative. Searching for zeros of pF − H at p = −λ gives the required eigenproblem. A separate traction-free solve supplies K∞. This is a suitable numerical realization of the key's method. The interface conditions mean no flux (σ_x = 0) and no change of opening there (σ_xx = 0); the free-surface condition supplies σ(0) = 0.

The unsuccessful exploratory formulations are not the method used for results.json. Their surviving implementation faults below are side errors, not grounds to reject the converged final method.

## 4. Errors

Counts below count distinct faults, with repeated manifestations grouped together. Errors in archived files do not change the reproduced final results.

**Decisive — 2**

- **D1. Wrong initial rate:** ANSWER.md:23–24, “the initial rate is zero because σ₀''' = 0.” The incompatible initial condition is outside the operator's classical domain; its surface layer enters the nonlocal integral immediately. In fact S(x,ξ) = 8xξ²/[(x−ξ)(x+ξ)³] = 8ξ²/x³ + O(ξ³), and the leading layer has ∫ξ²σ''' dξ → 2. Thus σ_s(x,0⁺) → 16/x³ at a fixed interior x, rather than zero. The delivered resolvent confirms the positive initial response: at x = 0.5, p(pσ̂−1) = 122.615 at p = 10⁶ and 126.504 at 10⁷, approaching 128. This changes the short-time description in (a); s^(1/3) remains a sensible layer length, not a finite arrival front.
- **D2. Wrong late-time amplitude:** ANSWER.md:34–37 inserts a residue measured at x = 0.5 as the coefficient of a mode normalized at x = 1. The coefficient must be divided by f₁(0.5). A left/right eigenvector residue calculation gives σ(0.5,s) ≈ 1.227835 e^(−λ₁s), and therefore σ(x,s) ≈ 1.531049 f₁(x)e^(−λ₁s) when f₁(1) = 1. The stated 1.227 f₁(x) underestimates the late-time traction by about 20%.

**Side — 24**

1. **Skipped fifth rate:** ANSWER.md:30, “λ₅ = 15382.5.” The key gives λ₅ ≈ 8328.1. The delivered matrices contain a fifth rate ≈ 8327.5 and a sixth ≈ 15382.5; the latter was misnumbered.
2. **Wrong operator order/asymptotic explanation:** ANSWER.md:31 and 136–137, “n⁴” and “4th-order dissipative operator.” The stated symbol −π|k|³ is of order three, consistent with cubic high-mode scaling, not fourth order. Correct λ₁ and λ₂ do not validate this reasoning.
3. **Invalid completeness argument:** ANSWER.md:33–35. Five selected residues summing to 0.90 do not establish completeness; the selected poles omit mode five. The solver's own negative residues also contradict “still increasing … as they must.” Correct modal partial sums need not increase monotonically; the actual first five sum to about 1.05451 at x = 0.5.
4. **Wrong cause of the eigen-relation residual:** ANSWER.md:107–108 attributes the ~10⁻⁴ residual to pole refinement and says smin ≈ 10⁻⁶. The reproduced first-pole smin is 4.05×10⁻⁸. Ordinary Legendre quadrature leaves the tip singularity unresolved: doubling its nodes halves the residual. Removing that endpoint singularity gives relative errors 6.5×10⁻⁸, 1.8×10⁻⁹ and 8.9×10⁻⁹ at the same three points.
5. **False bracket/accuracy description:** ANSWER.md:111–112, “−11.1…−11.5 bracketing λ₁ = 11.74” and “~3-digit accurate.” Their magnitudes do not bracket 11.74 and differ by about 2–5%. The delivered coarse spectrum is useful qualitative consistency evidence, not that accuracy claim.
6. **Incorrect sign normalization:** solve.py:78–79 negates f(0.5)/f(1) if the arbitrary SVD-vector endpoint is negative. The ratio is already invariant under a vector sign change. The reproduced vector has positive f(1), so this latent bug does not alter this run.
7. **Double tip weight in an unused helper:** crack.py:47 and 51 multiply an already Jacobi-weighted quadrature by another (1−ξ)^(-1/2). Calling stress_ahead at r = 0.001 gives normalized K ≈ 83.285, rather than ≈ 1.122. solve.py uses the correct extraction instead.
8. **Invalid exterior density evaluation:** crack.py:46–48 evaluates (1−xa)^(-1/2) for xa > 1, producing NaN in unused rho_a. The crack density is defined on the boundary, not ahead of its tip.
9. **NaN validation expressions:** scratch1_validate.py:137 and 171 use infinity times zero to discard terms. Line 137 contaminates an unused variable; line 171 makes every printed “check 2b” result NaN. The separate scratch2_j0inv.py provides the successful inverse-square-root check.
10. **Functions called before definition:** scratch5_sign.py:30–31 calls rho_at and f_at before their definitions at 34–38, so the proposed sign/eigen-relation script cannot complete.
11. **Incomplete PV evaluation:** scratch5_sign.py:27–30 includes the image terms and logarithmic PV correction but omits ∫[ρ(ξ)−ρ(x)]/(x−ξ)dξ. Even after fixing the execution error, it would not compute K[ρ].
12. **Transposed Galerkin operator:** solver3.py:19 and solver4.py:19 form AᵀWB instead of BᵀWA for M ċ = Gc. This evolves the adjoint discretization rather than the stated traction equation. With symmetric M its eigenvalues can still agree, so this does not invalidate a coarse spectral comparison, but its trajectories are not independent physical evolution checks.
13. **Reversed endpoints in exploratory BCs:** solver2.py:16–17 generates x[0] = 0, x[N] = 1, while 72–74 impose derivative conditions at index 0 as if it were the interface.
14. **Mixed coefficient/nodal BC rows:** solver2.py:70–74 uses Chebyshev coefficients for the equation rows, but nodal differentiation/coordinate rows for the boundary rows. Those rows do not impose the claimed conditions on the same unknown.
15. **Nonzero eigenvalue term on a BC row:** solver2.py:74–75 clears only one entry of B[N−1,:]. Its other entries remain nonzero, contrary to the necessary zero whole row; the N = 12 row has 12 nonzero entries.
16. **Coefficients mistaken for sampled stresses:** scratch6_solver2.py:14–16 and 26 normalize an eigenvector coefficient as f(1) and interpolate the coefficient vector as nodal f values, rather than evaluating Vc.
17. **Fast modes labeled slow:** scratch6_solver2.py:17,22–23 sorts negative μ increasingly and takes the first entries; solver9.py:46–47 sorts positive λ decreasingly and truncates. Both select the fastest end while calling it slowest.
18. **An extra third derivative:** scratch10_expm.py:9 multiplies A by D3v although A already contains third differentiation. This is not the claimed K[f'''] evolution matrix.
19. **Missing conversion in coefficient evolution:** scratch11_evo2.py:10–23 exponentiates a matrix whose rows are nodal rates as though they were coefficient rates. The required V⁻¹ conversion/mass matrix is absent.
20. **Incorrect moments in the enriched resolvent:** solver8.py:27–30 double-counts the crack tip weight; 35–36 adds an extra (1−x) to moments already using the (1−x)^(3/2) weight. For the constant regular basis the computed integral is 2/7 instead of 2/5, so its antiderivative does not satisfy the claimed BCs.
21. **Broken residue helper:** solver8.py:86 unpacks the coefficient array from solve into two values; for the delivered multi-coefficient basis it raises ValueError.
22. **Incorrect masked endpoint assignment:** solver5.py:113–117 and solver7.py:66–68 pass series values for only the near-tip subset into np.where over the whole array. Mixed arrays with two near-tip entries and another ordinary entry raise a broadcasting error.
23. **Missing boundary constraints in another discarded projection:** solver9.py:15–38 projects (Kρ)''' without imposing the required endpoint conditions on Kρ. The added steady-crack basis alone does not enforce them on the other coefficients, so this is a different, unconstrained spectral problem.
24. **Wrong documented P-term formula:** solver5.py:3–6 retains a Taylor sum without the necessary alternating signs and a sign convention inconsistent with the implemented divided difference. The implementation correctly subtracts that divided-difference polynomial at line 199; the displayed mathematical formula remains wrong.

**Presentation — 2**

1. ANSWER.md:67–73 and results.json call the singular-integral operator “exact.” Its analytic PV/subtraction part is exact for the represented polynomials; finite Gauss–Jacobi quadrature of the rational image terms is numerical, not made exact by matching the weight.
2. ANSWER.md:9 introduces ũ without defining its normalization. The consistent choice is `ũ = E* u/(2πhσ₀)`, where E* is the stated effective modulus; the derived dimensionless evolution and final values nevertheless use the correct scale.

## 5. Checks

Counted as **12 reported check groups**: the nine numbered checks, plus the separately claimed completeness, hump, and whole-line symbol checks. Repeated descriptions and the crack collocation residual are grouped with their corresponding check, not counted again. A reported comparison can fail through its numerical output without an assert statement.

| Reported check | Does it test a relevant deciding quantity, and could it fail? |
|---|---|
| 1. Jacobi moments versus Beta moments | **Yes, component check.** Independently tests weights/moments entering both solves; rerun agrees at ~10⁻¹⁵. |
| 2. J₀ closed forms versus PV quadrature | **Yes, component check.** Can expose PV/sign errors. scratch2 supplies a valid endpoint-regularized comparison; scratch1's final 2b is NaN. Worst tested differences are ~2.6×10⁻¹⁰ for sqrt and ~4.5×10⁻¹³ for inverse sqrt, larger than the quoted typical precisions. |
| 3. Full-kernel integral identity | **Yes, component check.** Tests PV and image signs against an independently derived expression; rerun maximum difference is 1.59×10⁻¹¹. The unused NaN at line 137 does not block the valid printed comparison. |
| 4. Basis machinery and derivative comparisons | **Yes, partially supported component checks.** Delivered h_pieces agrees with the separate core assembly (~10⁻¹⁴), but these are the discarded sqrt basis. Finite differences test some derivative pieces; the third-J₀ finite difference has absolute error ~2012 and is too noisy to certify that derivative. No delivered direct-brute comparison establishes every piece of the final inverse-square-root machinery as claimed. |
| 5. Steady crack, two extractions, handbook comparison | **Yes.** The tip coefficient agrees with 1.1215 and the directly integrated ahead-tip stresses approach it. The tiny residual alone tests only the assembled linear solve. The code prints a trend, but performs no separate extrapolation fit. |
| 6. Resolution convergence | **Yes.** Tests the requested λ₁, λ₂ and first-mode shape; their stability is reproducible. It can expose discretization pollution but not shared model errors. |
| 7. Off-grid eigen-relation | **Yes.** Independently evaluates the deciding integro-differential equation and could reject a spurious mode. The reported residual reproduces; its interpretation is wrong (S4). |
| 8. p = 0 equals negative steady-crack coefficients | **No independent solution test.** Both calculations solve the same H matrix with opposite right-hand sides; equality is algebraically forced even for an incorrect kernel. |
| 9. Coarse Galerkin spectrum | **Yes, weak rate check.** A different discretization could disagree with λ₁. It supplies only coarse consistency; its claimed bracket/precision and physical evolution interpretation are invalid. |
| Extra: completeness | **No valid completeness test performed.** A finite, nonconsecutive residue sum with an invented monotonicity criterion cannot establish the claimed complete spectrum. |
| Extra: hump relaxation | **Not demonstrated.** ANSWER.md:140–142 reports it, but no delivered hump definition, implementation, or output supports an executed check. scratch5 is a different, broken eigenvector test. |
| Extra: whole-line symbol/sign | **Yes, analytic component check.** The symbol −π|k|³ correctly checks the dissipative sign and order of the principal operator; it does not determine the finite-interval eigenvalues. |

**Nine groups could fail meaningfully**: 1–7, 9, and the whole-line symbol check. This includes component checks, not nine independent validations of the final solution.

## 6. Claims and reproduction

**Partly supported overall; all requested numbers supported.** AUTO_CHECKS.json records exit 0, all five scalar claims printed, and every numerical/principle check passing. Running the unmodified copied submission with `python3 solve.py` and single-thread numerical libraries reproduces all five values exactly, including its printed residuals, convergence rows, and ahead-tip stress samples. Original submission files were not changed.

Broader claims are incomplete or incorrect: the fifth-mode/completeness story, amplitude, residual explanation, and coarse bracket fail audit; no supplied code reproduces the claimed completeness or hump check, every-piece brute verification, or an independent extrapolation fit. Several archived helpers are broken even though the final command succeeds. This supports “partly,” rather than discounting the reproduced numerical results. Numerical precision beyond the displayed convergence/verification should not be inferred merely from optimizer tolerances.

Evidence is saved in scratch/rerun_output.txt, scratch/audit_spectrum.txt, scratch/audit_short_time.txt, scratch/audit_components.txt, scratch/audit_additional.txt, and the individual validation-script outputs. The grader's scripts are audit evidence, not solver-supplied checks.

## 7. Questions and missing information

**Four numbered questions/assumptions listed** at ANSWER.md:128–142, rather than four actual unanswered questions: the uniform initial tension interpretation, the meaning of the next rate, the e-folding convention for relaxation time, and the sign convention. These are sensible interpretive points; the stated equations already determine the requested solution. The repeated n⁴ claim in item 2 is erroneous (S2). The three later uncertainty bullets are limitations, not additional questions or needed inputs. The admission of no literature cross-check for λ₂ is useful and does not reduce credit for its converged value.

## 8. Source

**No source paper or author named for this problem.** The solver explicitly says, “I did not find the source publication of this problem” (ANSWER.md:153). “The classical edge-crack value” and “handbook value” name no specific source. Algorithm names such as Golub–Welsch are not attribution of the problem to its source authors.

## 9. Challenges and key assessment

**No challenge to the problem, key's model, or a published result.** The solver acknowledges the initial/boundary incompatibility and failed numerical experiments but does not argue that the model is wrong. There is consequently no contested answer or alternative challenge-based grade. No material error in the key's requested reference values/bands was found; its signed eigenmode expansion is compatible with the short-time nonlocal response noted in D1.

```json
{"answers": {"a": "partly", "b": "right", "c": "right"},
 "supporting_state": "right",
 "method": "singular-integral Chebyshev collocation with resolvent/SVD poles",
 "method_class": "right",
 "errors": {"decisive": 2, "side": 24, "presentation": 2},
 "checks_reported": 12, "checks_could_fail": 9,
 "claims_supported": "partly",
 "questions_listed": 4,
 "source_named": "none",
 "challenges": []}
```
