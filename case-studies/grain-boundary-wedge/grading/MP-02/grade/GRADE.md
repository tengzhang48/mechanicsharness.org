# Grading profile

## 1. Answers

| Question | Judgment | Reason |
|---|---|---|
| (a) | **Partly right** | Both requested rates, λ₁ = 11.7426680713 and λ₂ = 264.6342242, meet the key's exact-solution bands (1% and 2%). Surface relaxation and eventual first-mode dominance are correctly described, but the stated late-time traction formula omits the initial-condition-dependent modal amplitude. |
| (b) | **Right** | The limiting boundary is a traction-free, crack-like diffusion wedge, elastically equivalent to an edge crack. K∞/(σ₀√(πh)) = 1.1215222555 is inside 1.09–1.15. The extra finite-time statements have errors listed below. |
| (c) | **Right** | The reported asymptotic scale 1/λ₁ = 0.0851595220 t₀ is inside 0.05–0.15 t₀ and consistent with its own λ₁. The extra 1% threshold is inaccurate. |

All five required numerical outputs pass the key's bands. These judgments credit the numerical results while also assessing the written description.

## 2. Supporting state

**Right.** The requested normalized intermediate profile is f₁(0.5)/f₁(1) = 0.8019565030, within 2% of the key's 0.8020. The code uses the stated kernel and time unit, derives the correct evolution sign and eigenproblem λf = −L[f‴], and enforces f(0) = f′(1) = f″(1) = 0. The steady-density-to-K normalization is also correct. The copied simulation reproduces the normalized profile. The missing modal amplitude affects the absolute traction formula, not this normalized quantity.

## 3. Method

**Right method: singular-integral spectral collocation.** The solver uses a Chebyshev basis containing the tip power (1−Z)^(5/2), smooths the density singularity using Z′ = 1−y², subtracts the Cauchy pole analytically, and solves Bc = λAc. This is an appropriate implementation of the key's singular-integral route. A separate singular-basis steady solve produces the edge-crack limit. Polynomial collocation and matrix-exponential evolution provide additional comparisons.

The endpoint row uses the finite-part/interior limit of the transformed integral; an ordinary endpoint principal-value integral would diverge. The submission explicitly recognizes the finite-part interpretation in its kernel checks, and the resulting eigenvalues agree with the key. No change to the physical model is needed.

## 4. Errors

Count distinct errors once, including notation and baseline omissions. Counts: **1 decisive, 4 side, 3 presentation**.

1. **Decisive — missing first-mode amplitude.** `submission/ANSWER.md:80` states “σ_gb(z,t) ≈ σ₀ · f₁(z/h) · exp(−λ₁ t/t₀)” while fixing f₁(1)=1. For the uniform initial state, its own evolution gives a coefficient c₁ ≈ **1.531049**, not 1. At τ=0.1 the tip traction is approximately 0.473164 σ₀, whereas that formula gives 0.30905 σ₀. The correct late-time expression includes c₁. This changes the quantitative evolution in (a), although the rates and normalized shape remain right. See `scratch/audit.stdout.txt` and the basis-refinement results in `scratch/growth_convergence.stdout.txt`.

2. **Side — pointwise monotone growth asserted without being true.** `submission/ANSWER.md:92`: “The wedge grows monotonically”. The stated linear model has an early local decrease of opening. Defining W = E*u/(2πσ₀h), conservation gives ∂τW = −Σ″. At Z=0.889, τ=0.0025, the copied solver gives ∂τW ≈ **−1.95336**, consistently for N=24, 56, and 72. Thus the opening profile does not grow monotonically at every depth; this is not a claim about total inserted volume. The limiting wedge and K∞ are unaffected. See `scratch/growth_convergence.py` and its output.

3. **Side — crack-tip singularity incorrectly postponed to infinite time.** `submission/ANSWER.md:100`–101: “the √-singular crack field is attained only as τ→∞.” Bounded traction on the boundary does not rule out a stress-intensity field around its tip. The mode's Σ‴ ∝ (1−Z)^(-1/2), integrated through the conservation law, already creates a square-root opening and nonzero K at finite time. Using its own matrices gives **K(0.1t₀)/(σ₀√(πh)) ≈ 0.71029**. It is the traction-free limit and limiting K that are approached asymptotically. See `scratch/audit.py` and its output.

4. **Side — inaccurate 1% tip-traction threshold.** `submission/ANSWER.md:107` says “the tip traction is below 1% of σ₀ by τ ≈ 0.4”. Its delivered `sim.py` instead prints **Σ(1,0.4)=0.013966**, or 1.397%. A finer time grid crosses 0.01 at approximately **τ=0.4285**. This does not change the accepted e-folding scale in (c).

5. **Side — quadrature accuracy overstated.** `submission/ANSWER.md:50` claims “machine-precision values of L[g]”; lines 124–125 suggest differences much greater than 10^-10 would flag failure. Its own constant-function test has error **−1.38×10^-7** at Z=0.001. Even with the advertised M=384, testing at the first N=72 collocation point, Z≈0.00046294, gives constant-function error **−9.14×10^-8**. These errors are small enough for the grading bands, but are not machine precision and exceed that universal diagnostic threshold. See `scratch/additional_audit.stdout.txt`.

6. **Presentation — operational mean-traction baseline unstated.** `submission/ANSWER.md:105`–107 gives mean-halving and 1/e times near 0.085 and 0.11. `submission/sim.py:53` actually uses the mean at **τ=0.025** as the reference. With the initial mean as reference, the times are approximately **0.0655** and **0.0915**. The reported early-reference times are reproducible, but the reference time should be identified. This is a baseline omission, not an error in 1/λ₁.

7. **Presentation — transformed function changes name without definition.** `submission/ANSWER.md:46` defines G(y), but line 48 writes the subtraction formula using H(y). The expression is correct with H=G; that equality is unstated.

8. **Presentation — dimensionless opening undefined.** `submission/ANSWER.md:69` uses “Ã = (E*/2πσ₀) ∂_Z U” without defining U. The dimensions and the equality in `submission/steady.py:8` require **U=u/h**. The implemented K factor is correct.

The polynomial extrapolation claim discussed in item 6 below lacks supporting evidence; it is not counted here as an established false numerical value. The initial/boundary compatibility observation is also not an error.

## 5. Checks

Count the **six numbered check groups** in `submission/ANSWER.md:115`–170; subcases within a group are not counted separately. **All six could fail**, with the limits below. Printed diagnostics count as checks even though the scripts generally do not enforce assertions.

| Reported check | What it tests; ability to fail |
|---|---|
| 1. Kernel quadrature | Tests the deciding operator using an analytic constant-density integral, adaptive reference integration, and polynomial density. Wrong pole subtraction or image terms can cause discrepancies. It is useful, but its stated universal accuracy threshold is not met. |
| 2. Fractional-basis convergence | Directly tests stability of λ₁, λ₂, and f₁(0.5). Drift or oscillation could fail it. The entire quoted N=40–72 table reproduces. Convergence alone cannot detect a consistently wrong operator. |
| 3. Polynomial discretization | Tests λ₁ and λ₂ with a different basis and collocation grid. Disagreement can reveal basis or assembly errors. The quoted table reproduces and approaches the main values, but the precise extrapolated limits are not delivered. Both discretizations share `LOperator`, so this does not independently validate the kernel. |
| 4. Full spectrum and residual | Tests whether a slower complex/other mode was discarded and whether the eigenpair satisfies the assembled problem. These conditions could fail. Reconstructing N=30, M=384 confirms all 31 eigenvalues finite, the stated first eight real modes, complex-tail magnitudes ≥ approximately 4.49×10⁵, and relative mode-1 residual **8.16×10^-12**. A small residual alone only verifies the discrete equation. |
| 5. Time evolution | Tests the observed late decay and normalized profile, at N=40, against the main N=72 results. It could reveal a sign, evolution-assembly, or mode-resolution inconsistency. Both quoted outputs reproduce. This reuses the same elastic matrices and eigendecomposes the evolution generator; it is not an independent validation of the physical operator and does not test the absolute amplitude missing from the answer. |
| 6. Steady convergence/classical K | Directly tests K∞ and its normalization against a reference and increasing Nc. Wrong scaling or kernel can shift the result. The Nc=32,48,64,96 values reproduce. The key independently supports the approximately 1.1215 reference; the submission supplies no precise citation establishing its reference to seven decimal places. |

## 6. Claims and reproduction

**Partly supported overall; all five required numerical claims are supported.** `AUTO_CHECKS.json` records successful reproduction, 5/5 claims printed, and all key/principle checks passing. Rerunning `python3 solve.py` in `scratch/submission` exits successfully and regenerates all five values **identically**, including unrounded JSON values. `test_kernel.py`, `eig_poly.py`, `steady.py`, and `sim.py` also exit successfully, with outputs saved alongside the copied files.

The grader reconstructed the otherwise unprinted high-N convergence, full-spectrum/residual, and high-Nc steady checks; they support the quoted results. There is no mismatch in the required numerical outputs or the quoted check tables. There are, however, the narrative discrepancies in item 4, particularly the missing amplitude and the 1% statement contradicted by the simulation.

The statement at `submission/ANSWER.md:141` that “power-law extrapolation gives λ₁ → 11.7420 and λ₂ → 264.63” is **not demonstrated by the delivered code/output**: no extrapolation model, selected fitting points, fitting implementation, or fit output is supplied. `eig_poly.py` only prints finite-N eigenvalues. Fits to those values depend on the model and selected points; for example, a free three-parameter power-law fit to all printed points gives different limits. Credit the reproducible approach toward the main values, not the claimed precise extrapolation. The high-precision external fracture reference is likewise named but not precisely sourced.

## 7. Questions and missing information

**Zero questions or missing-information requests.** The “Questions / assumptions” section contains five assumptions, interpretations, or uncertainty notes, rather than questions. Taking the supplied kernel literally, using f₁(1)=1, defining relaxation as the asymptotic e-folding time, and treating the initial surface jump with the boundary condition enforced for t>0 are sensible. No missing physical input prevented the requested calculation.

## 8. Source naming

The solver names fracture-reference authors: **“Wigglesworth/Tada edge-crack constant 1.1215”** (`submission/ANSWER.md:98`). It does **not** name Gao, Zhang, Nix, Thompson, and Arzt, or the source paper on grain-boundary diffusion wedges identified by the key. Thus authors are named for the limiting fracture comparison, but the original diffusion-wedge source is not identified.

## 9. Challenges and the key

There is one **checkable compatibility qualification**, rather than a challenge to the numerical solution: “The initial uniform traction is incompatible with Σ(0,t)=0 at t=0 (a jump)” (`submission/ANSWER.md:180`). Its argument is valid if both conditions are required pointwise at the initial corner: σ₀>0 and zero cannot both hold there. The proposed interpretation—uniform initial traction in the interior and the surface condition for t>0—is the usual transient interpretation. It is checkable directly at the corner and by examining the projected solution as τ approaches zero from above.

This qualification is not counted as an error. It changes none of the required eigenvalues, normalized eigenmode, steady K, or e-folding scale: against the key and with this qualification, the judgments remain **(a) partly, (b) right, (c) right**. No requested answer is contested. There is no challenge to a published numerical result and no evidence that the key's answers or grading bands are wrong.

```json
{"answers": {"a": "partly", "b": "right", "c": "right"},
 "supporting_state": "right",
 "method": "Singular-integral spectral collocation with a tip-adapted Chebyshev basis",
 "method_class": "right",
 "errors": {"decisive": 1, "side": 4, "presentation": 3},
 "checks_reported": 6, "checks_could_fail": 6,
 "claims_supported": "partly",
 "questions_listed": 0,
 "source_named": "Wigglesworth/Tada edge-crack constant 1.1215",
 "challenges": [{"claim": "Uniform initial traction conflicts with zero surface traction at t=0; enforce the surface condition for t>0.", "checkable": true}]}
```
