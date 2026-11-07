# Submission grade

## 1. Answers

| Question | Judgment | Reason |
|---|---|---|
| (a) | **Partly right** | The eigenmode expansion, surface-first relaxation, and late-time dominance are correct, but λ₁=88.2491 and λ₂=1361.871 are wrong: the key gives 11.74265 and 264.6348, respectively; neither meets its accuracy band. |
| (b) | **Right** | The limiting plated wedge is elastically equivalent to a traction-free edge crack of depth h, and K∞/(σ₀√(πh))=1.1215 lies within 1.09–1.15. The solver explicitly identifies this as a handbook value rather than a computed elastic result. |
| (c) | **Wrong** | The definition t_relax=t₀/λ₁ is correct, but 0.01133 t₀ is outside 0.05–0.15 t₀; the correct rate gives about 0.08516 t₀. |

## 2. Supporting state

**Wrong.** The requested intermediate mode shape is f₁(0.5)=0.9705 instead of 0.8020, an error of about 21%, outside the 2% band. The qualitative boundary conditions and normalization are correct, but the reported profile is not the physical first eigenmode. No transient amplitudes or time-dependent intensity factors were required.

## 3. Method

**Weighted-Chebyshev singular-integral collocation; right method, incorrectly implemented.** The diffusion/elasticity combination, time normalization, eigenproblem λf=−PV∫Sf‴, and f(0)=f′(1)=f″(1)=0 are correct. Analytic treatment of the Cauchy term, paneled quadrature for the remaining kernel, triple integration, and a generalized eigenvalue solve are an appropriate numerical route. The implementation uses different signs for odd basis functions in its operator and traction matrices.

For s=sin²θ, **2s−1=−cos(2θ)**. Consequently, `solve.py:63` integrates Tₖ(1−2s), while `solve.py:53` and `:56` evaluate Tₖ(2s−1). Its implemented traction basis therefore satisfies fₖ‴=(−1)ᵏgₖ, although the operator assumes fₖ‴=gₖ. This is not a harmless common normalization.

The grader's [audit](scratch/audit_basis.py) corrects only the odd columns of the traction matrix at N=64. It obtains λ₁=11.7426677985, λ₂=264.6342230860, f₁(0.5)=0.8019565043, and t_relax/t₀=0.0851595240. All agree with the key. Correctly differentiating the submitted mode instead gives max|R[f]/f−88.24905|≈377.157 over the same 19 residual-test points, rather than 1.07×10⁻⁷. See [audit output](scratch/audit_basis.out).

## 4. Errors

Counts distinguish root implementation errors from distinct erroneous supporting statements. Repetitions and the four incorrect required numbers caused by D1 are recorded together, rather than counted again at every occurrence. Abandoned exploratory versions are not treated as additional final-solution errors.

| ID | Class | Evidence and error |
|---|---|---|
| D1 | **Decisive** | `solve.py:62–63`, contrasted with `:53` and `:56`; `ANSWER.md:88`: “gₖ ds = 2 Tₖ(cos2θ) dθ.” The missing minus sign in the Chebyshev argument reverses every odd traction-basis column relative to the operator. It changes both requested rates, the half-depth shape, and the relaxation time (`ANSWER.md:9–13`). The matrices solve a different equation. |
| S1 | **Side** | `ANSWER.md:22`: “the third is λ₃ ≈ 8331.” The physical third rate is about 1357.44. This unrequested rate is another consequence of the wrong spectrum. |
| S2 | **Side** | `ANSWER.md:65–66`, `:74`, `:132–136`: regularity purportedly “selects” the inverse-square-root behavior at both ends, and coefficient decay supposedly proves the exact analytic endpoint form. The three boundary conditions do not establish that exponent or analyticity; a basis weight permits endpoint cancellations. Convergence of this mismatched operator cannot validate the physical endpoint behavior. The weight remains usable once the sign error is fixed. |
| S3 | **Side** | `ANSWER.md:96`: “Every R[fₖ](Z) was compared.” `solve.py:124–125` tests only k=0,1,2,3,5 at six points, not every assembled entry of the 64×64 matrices. The comparison tests −PV∫Sgₖ, without checking that gₖ is the implemented fₖ‴. |
| S4 | **Side** | `ANSWER.md:119–120`: “a genuine eigenfunction of the continuum operator.” The small residual uses the same inconsistent `fk`/`Rf` pair as the solve (`solve.py:130–132`). It does not certify the physical eigenfunction. The code also filters out negative/complex eigenvalues (`:114–115`); the audit finds 30 negative real eigenvalues, including about −316.36 and −3845.24, in this erroneous matrix pair. Filtering is not evidence that the assembled spectrum has the physical positivity claimed at `ANSWER.md:22`. |
| S5 | **Side** | `ANSWER.md:138–141`: “Method independence,” using “both R and M scaled by ½.” Multiplying both matrices by the same nonzero factor preserves every generalized eigenvalue algebraically. The earlier code retains the same odd-basis sign error (`cheb5.py:24`). Agreement under this scaling is not an independent verification; the small numerical shift reflects its different quadrature. |
| S6 | **Side** | `ANSWER.md:26`: “monotone increasing.” The submitted converged mode slightly overshoots its endpoint value: `validate2.out:11` reports f(0.9)/f(1)=1.000009, while f(1)/f(1)=1 by definition. Thus even the computed profile is not strictly monotone. |
| P1 | **Presentation** | `ANSWER.md:27–28` lists 15 shape values for the 19 positions Z=0.05,0.10,…,0.95. The last four values are omitted, leaving the coordinate/value list incomplete. The code prints all 19. |
| P2 | **Presentation** | `ANSWER.md:38`: “there is no interface.” Equal elastic constants remove elastic mismatch; the film/substrate interface remains the locked, non-diffusing end of the wedge. The answer elsewhere correctly retains this tip condition. |
| P3 | **Presentation** | `ANSWER.md:108`: “a uniform Gauss grid.” The direct check uses nonuniform Gauss–Legendre nodes (`solve.py:87`). Its algebraic endpoint-quadrature error is plausible, but “uniform” is inaccurate. |

**Totals: 1 decisive, 6 side, 3 presentation.**

## 5. Checks

There are **four explicitly numbered checks**. Subcomparisons within a numbered item are counted together. “Could fail” means a nontrivial diagnostic capable of exposing some relevant errors, not a guarantee that it catches this submission's systematic sign error. None has an automatic pass/fail assertion.

| Reported check | Deciding quantity tested? | Could fail? |
|---|---|---|
| 1. Operator versus direct and paneled PV integration | Partly: tests singular-integral evaluation on five proposed derivative basis functions. It does not test their connection to the implemented traction basis. Both printed errors, 1.56×10⁻² and 4.37×10⁻¹¹, reproduce. | **Yes:** inconsistent Hilbert identities or quadrature assembly can produce disagreement. The implementations share the stated kernel, so agreement cannot independently validate that kernel's physical sign or factor. |
| 2. Eigenmode residual at 19 points | Tests the assembled equation away from collocation points, but not the physical equation because f‴ is misidentified. The printed 1.07×10⁻⁷ reproduces. | **Yes:** an off-grid discretization error can make this residual large; low-resolution exploratory outputs demonstrate that sensitivity. It misses D1. |
| 3. Basis/resolution convergence and coefficient decay | Directly compares the computed rates as N changes. `solve3.out` and `solve3b.out` support the quoted convergence values. It demonstrates numerical stability of the wrong equation, not correct rates or the asserted exact endpoint form. | **Yes:** inadequate resolution can cause drift. Coefficient decay is only an auxiliary diagnostic. |
| 4. “Method independence” under common matrix scaling | Does not independently test the physical rate, derivative reconstruction, or governing equation. | **No:** common scaling must preserve the generalized eigenvalues, whether correct or wrong. |

**Checks reported: 4; checks that could fail: 3.** No independent calculation of K∞ is claimed. The automatic inequalities λ₂>λ₁ and consistency with 1/λ₁ pass, but are not additional solver-reported checks and do not establish agreement with the problem.

## 6. Claims and reproduction

**Partly supported overall; numerical reproduction is fully supported.** `AUTO_CHECKS.json` records exit 0 and five of five numeric claims printed. The grader rerun on a copy also exits 0 and matches all five `results.json` numeric values and the quoted operator/residual diagnostics; see [rerun output](scratch/reproduction/rerun.out). The auxiliary convergence outputs support their numerical table. The emitted method string is abbreviated relative to `results.json`, but describes the same intended method.

The unsupported claims are the stronger interpretations of those checks: continuum-eigenfunction certification, exhaustive operator coverage, physical endpoint validation, and method independence. Reproducibility does not repair D1. The reproduction command itself prints one N=64 solve and its diagnostics; the resolution table and alternative normalization are documented in auxiliary scripts/outputs. K∞=1.1215 is assigned as a constant (`solve.py:120`), consistently with the solver's explicit admission that no elastic computation was performed.

## 7. Questions and missing information

**Four numbered entries, Q1–Q4** (`ANSWER.md:145–159`), plus two missing-work/uncertainty bullets (`:163–166`). Q1 sensibly interprets zero flux and the locked tip, but the conditions are already supplied. Q2 sensibly identifies the geometry dependence of K∞; the half-space is already specified. Q3 reasonably distinguishes an e-folding time from the nominal t₀, although its rate is wrong. Q4 is a numerical caveat rather than a missing physical datum; convergence precision does not resolve the wrong operator. The two additional bullets honestly acknowledge the borrowed crack factor and limited higher-mode precision. No genuinely missing information prevents solving the stated problem.

## 8. Source naming

**None.** The submission names no source paper or author. Its nearest attribution is `ANSWER.md:163`: “The 1.1215 edge-crack factor is a standard handbook value”; no handbook is identified. The authors listed in `KEY.md` are not credited to the submission.

## 9. Challenges and the key

**No challenge to the problem, physical model, key, or published result is made.** The questions concern assumptions, definitions, and numerical uncertainty. Rejecting unsuccessful numerical bases is not a challenge to the physical model. No answer is contested. The grader's checkable sign correction reproduces the key's rates, shape, and time scale; there is no evidence here that the key is wrong.

```json
{"answers": {"a": "partly", "b": "right", "c": "wrong"},
 "supporting_state": "wrong",
 "method": "weighted-Chebyshev singular-integral collocation",
 "method_class": "right",
 "errors": {"decisive": 1, "side": 6, "presentation": 3},
 "checks_reported": 4, "checks_could_fail": 3,
 "claims_supported": "partly",
 "questions_listed": 4,
 "source_named": "none",
 "challenges": []}
```
