# Grading profile

## 1. Answers

The operative references are the P374 corrections at KEY.md:18–33: K1_0 = 0.807593 and V1_short_coeff = 2.440220. The older paper values retained later in the key are not the grading targets. Nonzero values have a 1% band; K2_0 has absolute tolerance 0.01.

- **(a) Partly right.** The initial and final volumes and final K are correct; K1_0 = 1 earns the key's explicit partial credit for the undrained outer limit, but the asserted time independence of K is wrong, and the small-time volume coefficient is wrong.
- **(b) Wrong.** All three requested values are outside the bands. The negative long-time signs are correct, but replacing the poroelastic response by an elastic face traction −α gives incorrect magnitudes and an incorrect nonzero initial K.
- **(c) Partly right.** Four of nine numerical values meet the bands; K1_0 receives the conditional partial credit above, and the other four failed values are wrong.

All entries below use the problem's requested normalization.

| Quantity | Submission | Operative key | Judgment and reason |
|---|---:|---:|---|
| K1_0 | 1 | 0.807593 | Partly: correct for the reverse order of limits, not the drained-tip limit used by the key. |
| K1_inf | 1 | 1 | Right: the drained Griffith intensity. |
| V1_0 | 1.884955592 | 1.884955592 | Right: π(1−ν_u). |
| V1_inf | 2.513274123 | 2.513274123 | Right: π(1−ν). |
| V1_short_coeff | 0.708981540 | 2.440220 | Wrong: incorrect face contribution and omitted tip contribution. |
| p_faces_0 | 0.749063670 | 0.749063670 | Right: the pre-diffusion undrained pressure. |
| K2_0 | −0.89 | 0 | Wrong: exceeds the absolute tolerance by a large margin. |
| K2_inf | −0.89 | −0.33375 | Wrong: requires −η, with η = α(1−2ν)/[2(1−ν)]. |
| V2_inf | −2.236813969 | −0.838805239 | Wrong: requires −απ(1−2ν)/2, not −απ(1−ν). |

Under the key's alternative reading, K1_0 alone is right as an undrained outer intensity; (a) and (c) remain partly right because the short-time coefficient is still wrong. The submission does not identify or justify this order of limits.

## 2. Supporting state

**Right for the requested intermediate quantity.** The formula p_faces_0 = 2(ν_u−ν)/[α(1−2ν)] and its numerical evaluation are correct (ANSWER.md:40–42; solve.py:66). This is the induced pressure in the undrained solid before diffusion, rather than the pressure maintained on the drained boundary. The supplementary FEM is not a correct supporting solution of the stated equations; its defects are listed below.

## 3. Method

**Elastic effective-traction reduction with transient FEM; general_or_wrong.** The elastic slit formulas correctly give the two loading-1 volumes and final K. The extension to poroelasticity is invalid: although σ′ = σ + αpI is a legitimate definition, equilibrium gives ∇·σ′ = α∇p. Treating the skeleton as an elastic crack loaded only by its effective face traction omits this distributed forcing. Likewise, p = 0 on the faces does not make the evolving pressure field or K time independent.

The √t diffusion scale is appropriate, but the alleged weighted-profile integration is not shown and gives an incorrect coefficient. The key's diffusion/reciprocity calculation, including the corrected tip contribution, is needed. The supplementary full transient computation is unnecessary for the closed-form endpoints and does not repair these arguments.

## 4. Errors

Counts below count each distinct defect once; repeated statements in ANSWER.md, solve.py, and PROCESS.md are not extra errors. Supplementary FEM defects are side errors because solve.py computes the reported answers directly from formulas, without using the FEM results.

**Decisive: 5.**

1. **Incorrect constant-K argument and initial loading-1 intensity.** “The total-stress field — hence K — is time independent” (ANSWER.md:26–31; solve.py:49–52) wrongly infers the field from the drained face condition alone. It misses the drained tip region and produces 1 rather than 0.807593 under the key's limit. Partial credit for the alternative outer limit does not validate this reasoning.
2. **Incorrect small-time coefficient.** “A = 2√π (ν_u−ν)” (ANSWER.md:35–38; solve.py:60) gives 0.708982. Even the key's face-only contribution is 1.354055; adding the tip contribution 1.086165 gives 2.440220. Neither the required reciprocity calculation nor the tip contribution is supplied.
3. **Incorrect loading-2 initial response.** “K2_0 = K2_∞ = −α” in normalized units (ANSWER.md:45–50; solve.py:76) assigns a finite initial elastic traction response. The key requires K2_0 = 0. The effective-traction argument omits pressure-gradient forcing in the bulk.
4. **Incorrect loading-2 final intensity.** The same face-only replacement gives K2_inf = −α instead of −η (ANSWER.md:50; solve.py:77). The magnitude is too large by a factor 8/3 at the stated parameters.
5. **Incorrect loading-2 final volume.** V2_inf = −π(1−ν)α (ANSWER.md:51; solve.py:78) instead of −απ(1−2ν)/2 also overestimates closure by 8/3.

**Side: 8.**

1. **Incorrect element gradients.** fem_poro.py:55–58 multiplies row reference gradients by inv(J).T rather than inv(J). A linear patch audit gives grad(x) = [0,1] on the first triangle and [1,−1] on the second, rather than [1,0]. This corrupts elasticity, coupling, and diffusion, beyond the acknowledged finite-domain/mesh limitations.
2. **Incorrect fluid-content coefficients.** fem_poro.py:73–75 implements ζ = su div(u) + bcoef p with su = 1.068 and bcoef = 0. The stated two-dimensional laws require su = α = 0.89 and bcoef = α²(1−2ν)(1−2ν_u)/[2G(ν_u−ν)] = 0.23763. Thus the claimed FEM “exactly as given” (ANSWER.md:10,20) does not use the required ν_u or storage law.
3. **Inconsistent backward-Euler history term.** fem_poro.py:152–154 uses su Mp M_lump⁻¹ Cp u_n, whereas its left-hand displacement coupling requires su Cp u_n. The extra projection changes mass balance; the elastic-field audit finds a 5.20% relative discrepancy in that displacement history vector.
4. **Unjoined nominal crack tips.** fem_poro.py:18,23–26 splits the endpoint nodes as well as the interior faces, without tying endpoint displacements. The elastic solve has openings 0.26477 and 0.26734 at x = −L and L, respectively, where the slit opening must vanish. This also compromises the near-tip COD estimate.
5. **Permeability included twice.** fem_poro.py:102 puts kap into Kp, and line 137 multiplies Kp by kap again. The diffusion term would scale as κ². This is dormant for the delivered κ = 1, but is an error in the stated general FEM implementation.
6. **Loading-1 check misreads its own K history.** “Plateaus … i.e. the total-stress SIF is time-independent” (ANSWER.md:88–90) is an invalid implication: K(nu) increases from 0.87947 to 1.43016, by 62.6%. Agreement of the final plateau with an elastic calibration checks a final limit, not constancy in time.
7. **Volume-ratio check is not FEM confirmation.** ANSWER.md:88–91 identifies the trajectory with the required undrained/drained endpoints and presents 0.75 as confirmation. No t = 0 volume is computed; the first sampled/final FEM ratio is 0.56760. The stated 0.75 is simply the ratio of the answer's analytic formulas and does not validate the transient implementation or its initial limit.
8. **Loading-2 check misreads its own K history.** “A nearly time-independent K … K2_0 ≈ K2_∞” (ANSWER.md:100–101) contradicts K(nu) = −0.15114 initially sampled and −0.38227 finally: its magnitude grows by 153%. These finite-time samples also establish no t → 0 limit.

**Presentation: 2.**

1. **FEM K units and resulting calibration comparison.** fem_poro.py:163–168 calculates dimensional K, without dividing by √(πL), but lines 174–177 and ANSWER.md:69–76 compare it with normalized K = 1. At L = 1 the correct dimensional elastic reference is √π. The displayed 1.42905 corresponds to normalized 0.80626, about 19.4% below the reference, rather than “+43%.” The primary results.json normalization is correct.
2. **Pressure-location wording.** “On the drained crack line” followed by a nonzero undrained pressure (ANSWER.md:40–42; solve.py:63–66) conflates the imposed p = 0 boundary with the pre-diffusion adjacent solid. The requested pressure value itself is correct.

## 5. Checks

Count the three named composite checks A, B, and C once each: **3 reported, 3 could fail.** These are numerical comparisons capable of producing contrary evidence, even though the solver mishandles that evidence and the FEM is defective. No check establishes the corrected short-time coefficient or resolves the required tip/initial limits.

| Reported check | Deciding quantity tested? | Could fail? / finding |
|---|---|---|
| A: drained elastic calibration (ANSWER.md:69–78) | Tests FEM calibration against elastic V and K, not the deciding poroelastic limits. | Yes; it explicitly fails. V is 34.3% too small, and K has the additional unit mismatch above. No convergence repair is shown. |
| B: loading-1 transient and volume ratio (ANSWER.md:80–91) | Samples V and COD-derived K, but neither extrapolates t → 0 nor extracts the √t coefficient. The analytic volume ratio is not independent. | Yes; the K history actually contradicts the proposed constancy. The final elastic plateau and monotone volume alone cannot validate K1_0 or the growth coefficient. |
| C: loading-2 transient (ANSWER.md:93–103) | Tests the closure sign and samples COD-derived K; it does not verify the proposed factors or the initial limit. | Yes; closure signs agree, but K constancy fails in its own output. A negative sign cannot distinguish −α from −η. |

There are no assertions or independent analytical checks in solve.py; it evaluates and prints the proposed formulas.

## 6. Claims and reproduction

**Partly supported.** AUTO_CHECKS.json reports a successful solve.py rerun, all nine values printed, four value checks passing, and both volume-ratio/long-time-Griffith principle checks passing. My rerun in scratch/reproduction exits successfully and regenerates results.json exactly, including its formula strings. Thus there is no mismatch between the delivered primary numbers and their reproduction.

The separate unmodified fem_poro.py rerun also exits successfully and reproduces fem_out.txt exactly, including the excerpts pasted into ANSWER.md. However, these outputs do not support the claims of time-independent K, the independent volume-ratio confirmation, or verification by the full stated Biot model. The numerical magnitudes are explicitly disclaimed by the solver and therefore provide no corroboration of the erroneous coefficients. Detailed reproducible diagnostics are retained in scratch/audit.py and scratch/audit.txt.

## 7. Questions and missing information

**Three listed items** appear under “Questions / assumptions / uncertainty” (ANSWER.md:105–115), rather than three unanswered requests for information.

1. Loading-2 sign versus a hydraulically loaded crack: the distinction and superposition observation are sensible, but the task already specifies the two separate face conditions, so no clarification is needed.
2. Uncertain loading-2/short-time magnitudes: sensible given the failed calibration, but the accompanying assurance that K is time independent is unsupported.
3. Two-dimensional stress trace: a sensible consistency point and the correct reading of the stated indices. The FEM nevertheless uses incompatible fluid-content coefficients.

## 8. Source

Yes. ANSWER.md:10–12 names the “standard Rice–Cleary (1976) ‘crack-line’ reduction”; solve.py:4 also says “Rice & Cleary 1976 conventions.” No source paper title or equation is supplied, and this attribution does not establish the validity of the face-only poroelastic reduction.

## 9. Challenges and the key

**No challenge to the problem, key model, or published result is made.** The hypothetical discussion of hydraulic loading is an interpretation note, not a checkable argument that the specified loading is wrong. The constant-K assertions are ordinary solution claims, contradicted by the submitted computation, rather than well-argued challenges. No answer is marked contested.

The supplied key already records the initial-limit ambiguity and the corrections to the paper's K1_0 and short-time coefficient. Those qualifications are applied above. This audit supplies no independently established further key error.

```json
{"answers": {"a": "partly", "b": "wrong", "c": "partly"},
 "supporting_state": "right",
 "method": "elastic effective-traction reduction with transient FEM",
 "method_class": "general_or_wrong",
 "errors": {"decisive": 5, "side": 8, "presentation": 2},
 "checks_reported": 3, "checks_could_fail": 3,
 "claims_supported": "partly",
 "questions_listed": 3,
 "source_named": "Rice–Cleary (1976)",
 "challenges": []}
```
