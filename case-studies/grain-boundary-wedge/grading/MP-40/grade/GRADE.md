# Grading profile — MP-40

## 1. Answers

- **(a) Partly right.** The eigenmode description and late-time dominance are correct. `lambda1 = 11.871` is 1.093% above 11.74265: it receives the key's approximate-rate credit, but misses its 1% exact-rate band. `lambda2 = 257.493` is 2.699% below 264.6348 and fails the 2% band; the claimed limit near 258 and uncertainty of ±3 do not cover the reference rate. The argument for pointwise monotonic relaxation is also insufficient.
- **(b) Right.** A traction-free, crack-like wedge extending to the interface is the correct limit; `K_inf_norm = 1.1215` lies in the required 1.09–1.15 band.
- **(c) Right.** `relaxation_time_t0 = 0.0842` agrees with the solver's own reciprocal rate, 1/11.871 ≈ 0.08424, and with the key's 0.05–0.15 band.

These classifications agree with the numerical flags in `AUTO_CHECKS.json`. There is no contested answer.

## 2. Supporting state

**Right for the requested intermediate quantity.** The normalized first-mode traction at half-depth is 0.8001 versus the key's 0.8020, an error of −0.237%, within 2%. The delivered N=20 eigenvector reproduces that value and the displayed profile (`submission/out40b.txt:4–11`). At the orders actually used in `solve.py`, the basis satisfies the three homogeneous boundary conditions. The dimensionless evolution equation, kernel identity, and sign of the eigenproblem are correctly derived. Defects in other verification calculations are listed below and do not invalidate this accepted supporting value.

## 3. Method

**Right method: singular-integral spectral collocation with principal-value subtraction, followed by a geometric-tail extrapolation.** Combining diffusion and elasticity gives

`∂f/∂τ = N[f]`, `N[g] = −λg`, with `g(0)=0`, `g′(1)=g″(1)=0`.

The solver represents these conditions with `3x−3x²+x³` and `x(1−x)³ P_j(2x−1)`, forms `Oc = −λVc`, and solves using mpmath. This is the singular-integral collocation route appropriate to the key, although its numerical implementation and convergence evidence are inadequate for the requested second rate. The subtraction formula and its logarithm sign in the final code are correct. The edge-crack limit is identified physically and its standard factor is inserted as a constant, rather than calculated from the kernel.

## 4. Errors

Counts below identify distinct remaining errors in the delivered answer or verification code. Corrected exploratory failures, such as the earlier PV-log sign and eigenvector-row mistake, are not counted again. Approximate-rate credit is retained for λ₁; its explicitly qualified extrapolation is not treated as an exact solution.

**Decisive: 1.**

1. **The second decay rate is insufficiently converged and incorrectly bounded.** `submission/results.json:3` reports 257.493; `submission/ANSWER.md:79` says “λ₂ increases toward ≈ 258,” and lines 137–138 give “258 ± 3.” The key requires 264.6348 within 2%, i.e. approximately 259.342–269.927. Both the delivered number and its proposed asymptotic bound miss it. This changes answer (a).

**Side: 8.**

1. **False quadrature exactness.** “All integrals are exact for polynomial f with Gauss–Legendre” (`submission/ANSWER.md:44–45`; also `submission/solve.py:9–10`) is false: the reflected-kernel term is rational, not polynomial. In the local audit, even `φ₀'''=6` has an operator error of about −1.9963×10⁻⁶ at the first N=20 node, x≈0.004657. Exact integration changes the N=20 rates to approximately 12.176043 and 257.592173. The assertion of exactness is wrong; these small shifts do not independently change the answer classifications.
2. **The supposed alternative-basis check imposes different boundary conditions.** `submission/explore26.py:10–13` mixes ascending coefficients with descending-order `np.polyadd`, constructing `φ₀=x−x²+x³`. It has `φ₀′(1)=2` and `φ₀″(1)=4`, whereas both must vanish. The 24.8/25.0 rates in `submission/out30.txt` therefore do not establish the intrinsic “basis-sensitive” behavior claimed in `submission/ANSWER.md:109–112`. This erroneous check is not used to generate the final numbers.
3. **Incorrect diagnosis of the large-N collapse.** The claim that `V` becomes catastrophically ill-conditioned “even at 130 digits” and therefore only N≤20 can be used (`submission/ANSWER.md:122–128`) overlooks the float64 construction of monomial coefficients before high-precision evaluation (`submission/solve.py:31–62`, shared with the convergence probes). At N=24 the rounded basis already violates the boundary conditions: maximum column defects in value, first derivative, and second derivative at x=1 are about 35.91, 1995.04, and 11217.48. Increasing mpmath precision cannot restore lost input coefficients. A local calculation using exact coefficients and analytically integrated kernel remains stable through N=40, giving λ₁≈11.96910 and λ₂≈260.79163. This refutes the claimed intrinsic barrier, without claiming that N=40 is a fully converged solution. The selected N=16/20 basis coefficients themselves are exact.
4. **The claimed verification of an entirely real spectrum is false.** “All real (verified)” (`submission/ANSWER.md:82`) is unsupported by code that discards nonreal eigenvalues (`submission/solve.py:103–105`). The local rerun of its N=20 matrix finds six nonreal eigenvalues, including −191108.628±15225.636i. This does not disprove the real continuum spectrum in the key; it disproves the asserted numerical verification.
5. **Positive real decay rates alone do not imply pointwise monotonic traction.** The inference “so σ_gb relaxes monotonically” (`submission/ANSWER.md:82–83`) requires additional information about modal coefficients or positivity of the evolution. For example, `e^(−τ)−0.5e^(−10τ)` initially increases despite having positive real decay rates. Late-time exponential decay is correctly established, but this reasoning does not prove the stronger transient claim.
6. **Contradictory closed-form arithmetic.** `submission/ANSWER.md:129–130` allows `π⁴/8` near the numerical estimate but excludes `π(π/2)³` as “>13.” They are identical and equal approximately 12.176136, also correcting the printed 12.182. This is an error in a discarded candidate discussion, not in the delivered numerical computation.
7. **The boundary conditions do not themselves impose an interface maximum.** “It makes σ(1) a free maximum” (`submission/ANSWER.md:146–147`; similarly lines 100–102) confuses a free value with a maximum. Zero first and second derivatives alone do not establish a maximum. The computed first-mode profile does have the stated maximum; its attribution solely to the boundary conditions is the error.
8. **Supplementary verification is not packaged as reproducible within the submission.** For example, `submission/explore33.py:1`, `explore37.py:1`, and `explore40.py:1` load absolute paths under `<round-runs>/runs/MP-40/`, outside the supplied folder. They need path repairs to run from the delivered copy. The standalone `solve.py` works, but “Reproduce every number” with that command (`submission/ANSWER.md:3–6`) overstates its output: it does not execute the supplementary checks or print the complete convergence/profile tables.

**Presentation: 4.**

1. “σ₀ then decays” (`submission/ANSWER.md:35`) should refer to σ_gb; σ₀ is the fixed initial stress scale.
2. “f″ = −λg” (`submission/ANSWER.md:55`) is a mistyped eigen-equation; the actual equation is `N[g]=−λg`. The adjacent matrix formulation and code use the correct operator.
3. The subtraction limit “−f'''(x)” in `submission/solve.py:9` should be `−f''''(x)`. The implemented difference quotient and the answer's `−p′(x)` are correct.
4. `τ_relax = 1/λ₁ ≈ 0.084 t₀` (`submission/ANSWER.md:91`) mixes dimensionless τ with dimensional time. Correctly, `τ_relax≈0.084`, or `t_relax≈0.084 t₀`; the JSON value has the intended units.

## 5. Checks

Seven distinct verification claims are counted, including explicit claims outside the “Checks” section. The duplicated PV comparison is counted once. “Could fail” means capable of detecting a relevant implementation or numerical error; it does not mean an independent validation of the complete continuum solution.

| Reported check | Does it test the deciding quantity? | Could fail? / assessment |
|---|---|---|
| Symbolic identity of the two kernels (`ANSWER.md:28`) | Tests the operator reduction, not the rates. | **Yes.** The identity is correct; algebraic verification can expose a wrong kernel. |
| Quadrature versus independent symbolic PV integrals (`ANSWER.md:46–49,106–108`) | Tests the PV operator on monomials x⁴…x⁸ at x=0.4,0.7. It does not test endpoint accuracy or high-degree assembly. | **Yes.** `out33.txt` documents failure before the log-sign correction; `out33b.txt` documents agreement afterward. The retained code and local audit support these limited checks. |
| Entirely real, ordered positive spectrum “verified” (`ANSWER.md:82`) | Claims a spectral property, but only selected real, decaying eigenvalues are retained. | **No.** Selection guarantees the reported subset passes; excluded complex values are not tested and are actually present. |
| Backward Euler and eigenvector-convention comparison (`ANSWER.md:109–112`) | Tests rates and eigenvectors of the separate discrete matrix, not the final physical eigenproblem. | **Yes, limited.** `explore30.py:35–62` can expose indexing or integration mistakes. Its time-step rate is properly converted by `(exp(rΔτ)−1)/Δτ`, but agreement shares the same operator and the wrong boundary basis. It does not confirm 11.871. |
| N-convergence study (`ANSWER.md:65–79,113–115`) | Directly tests the numerical decay rates. | **Yes.** The low-order table is reproduced in the logs. The larger study actually fails to stabilize; the four selected rows and hardcoded mean difference ratio are not a convergence certification. |
| Mode-shape reproduction (`ANSWER.md:116–118`) | Tests the requested normalized half-depth traction. | **Yes, limited.** `out40b.txt:7,11` gives 0.800103 and the earlier Chebyshev/constraint implementation gives 0.800149 (`out25b.txt:3,8`), both rounding to 0.8001. These are different implementations of essentially the same polynomial collocation space, not independent continuum methods. |
| Repeat at 130 digits (`ANSWER.md:123–125`) | Tests sensitivity to mpmath arithmetic precision, not accuracy of the float inputs, quadrature, or boundary conditions. | **Yes, limited.** It could reveal precision sensitivity. `out39.txt` supports the observed repeatability, but the inference that it rules out input rounding is wrong. |

**Seven checks reported; six could fail in the limited senses stated.** No successful check demonstrates convergence of λ₂ to the accuracy requested, and no direct numerical K-limit check is reported.

## 6. Claims

**Partly supported overall; all five principal numerical output claims are supported.** `AUTO_CHECKS.json` reports exit status 0 and 5/5 numeric claims printed. A fresh run of `python3 solve.py` on the copy in `scratch/reproduction/` also exits 0 and exactly reproduces 11.871, 257.493, 0.8001, 1.1215, and 0.0842 after the delivered rounding. Its captured output is in `scratch/reproduction-rerun.txt`. There is no mismatch between these outputs and `submission/results.json`.

The code derives λ₁ from the N=16/20 rates with the fixed ratio 0.6875, takes λ₂ and the mode shape from N=20, and sets K=1.1215 as an explicit edge-crack constant. Thus reproduction is not independent confirmation of physical correctness. Logs support the limited PV, mode-shape, time-stepping, and precision observations; they do not support the claimed all-real spectrum, intrinsic high-order barrier, or quadrature exactness. The normal reproduction command does not run those supplementary checks. Local audit evidence is retained in `scratch/audit_basis.py`, `scratch/audit_basis.out`, `scratch/audit_spectrum.py`, and `scratch/audit_spectrum.out`.

## 7. Questions

**One labeled question/assumption item; no request for additional information.** The “Questions” section states: “The problem fixes ∂σ/∂z = ∂²σ/∂z² = 0 at z = h (rather than σ(h)=0). I used it verbatim” (`submission/ANSWER.md:146–147`). Using the specified conditions is sensible and requires no clarification. The attached assertion about a maximum is addressed under errors. Section 5 acknowledges numerical uncertainty and incomplete convergence, rather than identifying missing problem data; these are not counted as additional questions.

## 8. Source

**None named.** The submission names no source paper or author. “The standard edge-crack value” (`submission/ANSWER.md:88`) and “the classic edge-crack factor” (line 140) are descriptions, not source attributions. The paper and authors named in the grading key are not credited to the solver.

## 9. Challenges

**None.** The submission does not argue that the problem, its physical model, or a published result is wrong. Its discussion of slow convergence, possible endpoint logarithms, and an alternative smooth-tip interpretation concerns its own numerical uncertainty. The smooth-tip sentence (`submission/ANSWER.md:140–142`) supplies neither a demonstrated contradictory result nor a runnable alternative, and endorses the usual crack limit. No answer needs dual grading. The local audits expose defects in the submitted numerical work; they provide no evidence that the key is wrong.

```json
{"answers": {"a": "partly", "b": "right", "c": "right"},
 "supporting_state": "right",
 "method": "PV singular-integral spectral collocation with geometric-tail extrapolation", "method_class": "right",
 "errors": {"decisive": 1, "side": 8, "presentation": 4},
 "checks_reported": 7, "checks_could_fail": 6,
 "claims_supported": "partly",
 "questions_listed": 1,
 "source_named": "none",
 "challenges": []}
```
