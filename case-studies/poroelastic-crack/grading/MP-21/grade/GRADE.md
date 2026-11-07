# Grading profile

Graded against the **P374 corrected key**, using 1% bands for nonzero quantities and 0.01 for `K2_0`. The historical paper values retained in KEY.md are not the grading targets. The submission was reproduced on a copy in `scratch/reproduction/`; independent review calculations are in `scratch/review_checks.py` and `scratch/review_checks_output.txt`.

## 1. Answers

| Question | Judgment | Reason |
|---|---|---|
| (a) | partly | Both volumes and the long-time K are right; initial K uses the reverse-order, undrained outer limit, and the short-time volume coefficient misses both part of the face response and the leading tip contribution. |
| (b) | partly | `K2_0 = 0` is right; both long-time zeros are wrong because pointwise uniform pressure does not imply zero stress or zero crack opening. The challenge is assessed in section 9. |
| (c) | partly | Five of nine values pass the stated bands; `K1_0` receives the key's explicit partial credit for the alternative order of limits. The remaining three failed quantities are wrong. |

| Quantity | Submitted | Corrected key | Judgment |
|---|---:|---:|---|
| `K1_0` | 1 | 0.807593 | partly: outside the band, but correct as an undrained outer intensity |
| `K1_inf` | 1 | 1 | right |
| `V1_0` | 1.8849555922 | 1.8849555922 | right |
| `V1_inf` | 2.5132741229 | 2.5132741229 | right |
| `V1_short_coeff` | 1.1283791671 | 2.440220 | wrong |
| `p_faces_0` | 0.7490636704 | 0.7490636704 | right |
| `K2_0` | 0 | 0 | right |
| `K2_inf` | 0 | −0.33375 | wrong |
| `V2_inf` | 0 | −0.8388052385 | wrong |

The formulas for the five correct values agree with the key. The incorrect numerical values come from the same incorrect formulas and assumptions, rather than evaluation mistakes. The corrected tip limit is taken at fixed positive time before taking time to zero; the solver neither identifies nor analyzes that limit.

## 2. Supporting state

**Right for the requested supporting quantity.** The induced pressure before face drainage is correctly obtained from ζ = 0 and σ_xx + σ_yy = −2:

`p_faces_0 = 2(ν_u−ν)/[α(1−2ν)] = 0.7490636704`.

The outer undrained opening and its integral are also correct. The local constrained-strain diffusion calculation correctly obtains diffusivity c and its one-dimensional pressure profile. Its conversion into the full crack-volume increment, and the unrequested loading-2 zero-stress state, are incorrect as detailed below.

## 3. Method

**Method:** elastic Griffith limits, local oedometer diffusion, and scalar diffusion limits. **Class:** general_or_wrong.

Elasticity with ν_u and ν gives the correct loading-1 volumes and drained long-time K. The local one-dimensional problem is a useful component calculation, but its displacement is not the complete crack response. The needed short-time route combines diffusion with reciprocity and includes the tip contribution; the needed initial K requires the drained inner tip limit. For loading 2, a scalar pressure calculation plus an assumed freely swelling state omits the mechanical response to the expanding pressurized region. The numerical checks do not supply those missing solutions.

## 4. Errors

Counts below distinguish independent reasoning errors; repeated formulas and downstream numerical manifestations are not counted again. The accepted pressure-limit objection in section 9 is not an error.

**Decisive: 4.**

1. **Wrong initial K limit.** `submission/ANSWER.md:9,33`: “ζ ≡ 0 everywhere” and “K = √(πL) (SIF of the Griffith crack depends only on face tractions).” This treats the undrained outer stress field as the fixed-positive-time tip field. It gives 1 rather than 0.807593. The outer interpretation earns partial credit, but the reasoning does not establish the requested tip limit. The claim “SIF independent of time” in `submission/results.json:13` repeats this error.
2. **Incomplete face-volume response.** `submission/ANSWER.md:35`: “ΔV = 4L·δu_y ⇒ A = 8(ν_u−ν)/(√π(1−ν)).” The local strain integral assumes the outer displacement increment is zero and is treated as the entire face displacement. Reciprocity instead gives the face contribution `16(1−ν_u)(ν_u−ν)/[√π(1−ν)] = 1.354055001`, already different from 1.128379167 before adding the tips. Independently, eliminating p gives `e = D_u σ + η_u ζ I`, with `η_u=(ν_u−ν)/[α(1−2ν)]`; applying reciprocity to the face leak-off gives `A_faces=32η_u²κG/(c√π)`, the same 1.354055001.
3. **Wrong order of the tip contribution.** `submission/ANSWER.md:35,89`: “Tip-region contributions are O(ct/L²), subleading.” The corrected key has a leading √(ct/L²) tip contribution of 1.086165. Omitting it changes the requested leading coefficient, independently of error 2.
4. **Wrong stress selection in loading 2.** `submission/ANSWER.md:39,78,90`: “With p uniform, σ ≡ 0” and “uniform swelling ⇒ K2_inf = V2_inf = 0.” A compatible freely swelling local state is not thereby the limit of the initially unpressurized infinite solid. The remote transition region leaves a confining stress and the crack needs a traction-cancelling elastic perturbation. This changes both endpoints to `K2_inf = −α(1−2ν)/[2(1−ν)]` and `V2_inf = −πα(1−2ν)/2` in the requested units. Section 9 gives checkable evidence rather than relying only on disagreement with the key.

**Side: 5.**

1. **Wrong complex square-root branch in the Westergaard check.** `submission/solve.py:79–89` uses the principal `sqrt(z*z−L*L)`, which does not approach z in the left half-plane. At `(−100,0.5)` it gives σ_yy ≈ −2.000050 rather than decay to zero; at `(100,0.5)` it gives +0.000050. Normal stresses lack the required reflection symmetry. The right-side samples and the separate dislocation calculation still support the correct reported elastic volumes and drained K.
2. **Wrong extra hydraulic-loading results.** `submission/ANSWER.md:28`: “K = P√(πL) at both limits” and the stated drained elastic long-time volume for loading 1 + loading 2. Using the corrected endpoints, the combined normalized K values are 0.807593 and 0.66625, and the long-time V is 1.6744688844, not 2.5132741229. These quantities were not asked for separately.
3. **Scalar heat evolution presented as poroelastic verification.** `submission/solve.py:329–369` evolves only `p_t = ∇²p`. Compatibility and mass balance instead give harmonic `h = σ_kk + α(1−2ν)p/(1−ν)` and `p_t = c∇²p − [(ν_u−ν)(1−ν)/(α(1−2ν)(1−ν_u))] h_t`. The check never solves for h. Its printed values are valid for its scalar model, not demonstrated transient Biot pressures. This does not invalidate the separate mathematical objection to a bounded decaying steady pressure field.
4. **Overstated coverage of check 2.** `submission/ANSWER.md:65` says the check could catch “a factor-2 error in the two-faces × 2L bookkeeping.” `submission/solve.py:250–262` compares one face's local displacement and pressure profile; it does not form a numerical crack-volume increment and compare it with the reported A. Such a bookkeeping error could leave all those comparisons unchanged.
5. **Wrong explanation of the DD volume error.** `submission/ANSWER.md:57,91` attributes the discrepancy to an “endpoint clamp.” The dominant error is the 801-point trapezoidal integration of the square-root opening (`submission/solve.py:176,187`). Applying that same trapezoid directly to the exact drained ellipse gives 2.51315653, close to the DD result 2.51315776; the exact value is 2.51327412. The small discrepancy is within tolerance and does not change an answer.

**Presentation: 1.**

1. `submission/ANSWER.md:65`: “Ratios < 1 at t ≤ 1e−5” contradicts `submission/run_output.txt:32–36`, where the ratios are 5.674836, 3.275294, 1.987365, 1.356340 and 1.107215, all **above** 1. Grid under-resolution is real, but the stated direction is reversed.

## 5. Checks

Count the solver's **three numbered checks** once each, treating 3(a) and 3(b) as components of check 3. Checks 1 and 2 could expose an error in a requested elastic formula or the local contribution used to derive A. Neither verifies every part of the solution. Check 3 has a falsifiable pressure comparison, but cannot expose an error in its claimed mechanical endpoints.

| Reported check | Does it test the deciding quantity? | Could it fail? |
|---|---|---|
| 1. Westergaard / distributed dislocations / analytic elastic crack | Yes for the two elastic volumes, drained K and outer face stress; no for the drained inner initial tip K. | Yes: independently computed dislocation volumes, stress extraction and kernel tests can disagree with the elastic formulas. The square-root branch defect escapes the chosen samples. |
| 2. Local drainage FD / erf / surface displacement | Tests the local profile and displacement contribution, not the complete A including global mechanical adjustment and tips. | Yes for an incorrect local profile or strain-integral formula, as the unresolved-grid ratios demonstrate. A wrong full A, tip omission or face-count factor can still pass. |
| 3. Laplace box rejection of g, and parabolic pressure probes | Tests a scalar pressure candidate and finite-time scalar diffusion, not K2_inf or V2_inf. | A pressure candidate can fail, but the claimed zero K2/V2 cannot fail through this check: neither is computed or compared. The finite-box transient cannot establish the infinite-time mechanical state or an infinite-domain pressure asymptote by itself. |

**Checks reported: 3; checks that could expose an error in a requested result or its computed contribution: 2.** Only check 1 directly verifies requested final quantities in a physically appropriate elastic limit. Check 2 could catch, for example, an incorrect erf integral used in A, but cannot catch the missing global and tip contributions. Check 3 is excluded from this count because its pressure tests leave both wrong mechanical zeros untouched. There are no automatic pass/fail assertions in the script; failure here means a numerical comparison contradicts its stated prediction.

## 6. Claims and reproduction

**Partly supported overall; all numerical artifacts reproduce.** AUTO_CHECKS.json reports exit 0 and all nine scalar claims printed. The copied rerun also exited 0 with empty stderr; its `results.json`, including formulas, is identical to the delivered file, and its complete output is byte-for-byte identical to `submission/run_output.txt`. The displayed numerical excerpts match, except for the prose inequality identified above.

The automatic value checks and both principle checks are confirmed: five values pass, four fail, the volume ratio is 4/3, and long-time loading-1 K is 1. Reproducibility does not validate the physical formulas. The claim that the FD layer establishes the full A, and especially that pressure probes “pin down” zero long-time stress/opening, is unsupported: no coupled transient stress, K2 or V2 is computed. Those endpoint zeros are assigned in `submission/solve.py:68–69`.

## 7. Questions and missing information

**Four numbered questions/assumptions** are listed in `submission/ANSWER.md:80–85`. The loading-2 distinction between a global steady state and a pointwise infinite-time limit is sensible and important. The total-stress convention follows the problem. The observation about negative opening and absent contact is appropriate for the stipulated linear crack model. Using the supplied two-dimensional constitutive laws without another reduction is appropriate.

The three uncertainty bullets at lines 89–91 are disclosures, not additional questions. They correctly admit the absence of full two-dimensional transient verification, although the claimed order of the neglected tip correction and the explanation of quadrature error are wrong. No missing material information prevents answering the problem. The important initial-K order-of-limits ambiguity is not listed.

## 8. Source

The solver names a constitutive attribution: **“Rice–Cleary plane-strain form”** (`submission/ANSWER.md:7`; also “Rice-Cleary form” in `submission/solve.py:5`). Griffith and Westergaard are named for classical elastic solutions. No source paper for this particular crack problem is identified, and neither Detournay nor Cheng is named. The JSON records the literal attribution, not a citation to the key's source paper.

## 9. Challenges and key review

**Challenge 1: no bounded decaying loading-2 steady pressure field in an infinite two-dimensional domain.** Quote (`submission/ANSWER.md:82`): “has no steady state with p → 0 at ∞: p → 1 pointwise.” The argument uses two-dimensional recurrence, increasing interior pressures with increasing Laplace-box size, and transient pressure probes; it also rejects the unbounded candidate `1−Re[z/√(z²−L²)]`. This is checkable by the delivered code. The bounded-steady-state objection is well argued and is **not counted as an error**. Finite-box calculations alone do not prove an infinite-domain asymptote, and the transient script is uncoupled, but the steady-state objection is also mathematically sound. It does not contradict the key's confining-stress result, which need not assume that rejected harmonic candidate.

**Challenge 2: the pointwise pressure limit selects free swelling and therefore zero long-time K and V.** Quote (`submission/ANSWER.md:39`): “zero relative displacement across the crack: K2_inf = 0, V2_inf = 0.” The argument is that p = 1, σ = 0 satisfies local constitutive equations, equilibrium and face traction, and that uniform swelling is compatible. This is checkable, but the delivered pressure checks do not test the required mechanical selection, and independent checks refute the inference:

- For an infinite solid with p = 1 inside a disk of radius R and p = 0 outside, displacement/traction continuity gives `u = a r` inside and `u = a R²/r` outside, with `a = α(1−2ν)/[4G(1−ν)]`. Every finite-R solution decays at infinity, yet at any fixed interior point, as R grows, pressure is 1 and stress remains `−α(1−2ν)/[2(1−ν)] = −0.33375`. The strain is 0.166875, not the claimed free-swelling strain 0.267. This exact counterexample shows why a distant transition region cannot be discarded when selecting stress.
- A check specific to crack volume uses drained elastic reciprocity: with zero mechanical face traction, `V = β ∫p S_aux dA`, where `β = α(1−2ν)/(2G)` and S_aux is the stress trace of the unit-pressure Griffith crack. In elliptic coordinates, its angular integral including the Jacobian is `−2πL² exp(−2μ)`. For pressure 1 inside μ < M and 0 outside, `V = −βπL²(1−exp(−2M))`. At M = 8 this is −0.8388051441, approaching the key's −0.8388052385 rather than zero. These are exact static/reciprocal checks of the claimed implication, not a new transient simulation.

The expanding pressurized region's confining stress, cancelled on the crack by an elastic face perturbation, also yields the key's `K2_inf = −0.33375`. Against the key, (b) and (c) are partly right. Against the accepted pressure-limit objection alone, they remain partly right. If the solver's additional free-swelling assumption were imposed, its two zero endpoints would be right by construction, making (b) right under that added model and leaving (c) partly right; that assumption is not established by the problem or the computation. Thus the **requested numerical challenge is not well argued mechanically and is not marked contested**. The sound pressure objection is preserved separately and unpenalized.

**Key review:** The P374 corrections govern grading; the older 0.75 and 1.354055 rows are explicitly historical. There is a separate dimensional transcription error in the historical opening–K relation at `KEY.md:61`: `√(πG)` should be `G√π`. From the drained elastic opening, `lim D_n/√r = 2(1−ν)K√(2/π)/G`, hence `K = G√π lim(D_n/√r)/[√8(1−ν)]`. The printed version has the wrong stress units. This typo changes neither the supplied normalized grading values nor this submission's judgments under either reading. No new evidence here overturns the corrected numerical key.

```json
{"answers": {"a": "partly", "b": "partly", "c": "partly"},
 "supporting_state": "right",
 "method": "elastic Griffith limits, local oedometer diffusion, scalar diffusion limits", "method_class": "general_or_wrong",
 "errors": {"decisive": 4, "side": 5, "presentation": 1},
 "checks_reported": 3, "checks_could_fail": 2,
 "claims_supported": "partly",
 "questions_listed": 4,
 "source_named": "Rice–Cleary plane-strain form",
 "challenges": [{"claim": "Loading 2 has no bounded decaying steady pressure field in infinite 2D; distinguish the pointwise limit", "checkable": true}, {"claim": "Pointwise uniform pressure implies free swelling and zero long-time K2 and V2", "checkable": true}]}
```
