# Submission grading profile

## 1. Answers

Use the **corrected P374 key** (`KEY.md:11–34`), not the superseded paper formulas retained later in that file. All nine reported numbers pass the stated bands: 1% for nonzero values and absolute tolerance 0.01 for `K2_0`.

| Question | Profile | Reason |
|---|---|---|
| (a) | **Contested; right against the corrected key and under the solver's challenge** | The undrained/drained volumes and final intensity are correct; the actual-tip spectral limit gives `K1_0 = 0.808`, and complete short-time growth gives `A = 2.44`, agreeing with the corrected `0.807593` and `2.440220`. General expressions specify the spectral system and order of limits (`submission/ANSWER.md:59–91`). |
| (b) | **Right** | `K2_0 = 0`, `K2_inf = -eta`, and `V2_inf = -pi*alpha*(1-2nu)/2` have the correct magnitudes and closure signs; the background-stress cancellation argument supports them (`submission/ANSWER.md:45`). |
| (c) | **Contested; right against the corrected key and under the solver's challenge** | Every required numerical value is within its band; the two challenged quantities differ from the corrected key by only 0.05040% and 0.009016%, respectively. |

The contested labels record the well-supported challenges discussed in section 9, not numerical failures against the current key. Against the **old paper-based reference**, (a) and (c) would be partly right: seven numbers agree, but the two disputed quantities do not. Against the **corrected key**, both are fully right. Under the solver's alternative initial-time interpretation, the outer undrained intensity is 1, while the reported actual-tip intensity remains 0.808; the solver explicitly distinguishes these readings (`submission/ANSWER.md:27,53`).

## 2. Supporting state

**Right.** The requested induced outer face pressure is `0.7490636704119851`, derived from the initial face stress trace `-2` and `zeta = 0`. The solver correctly distinguishes this undrained solid-side value from the imposed boundary pressure zero at positive time (`submission/ANSWER.md:43,50`). The initial opening uses `nu_u`, the final opening uses `nu`, and full opening rather than one-face displacement is integrated for volume.

Independent substitution verifies `eta = 0.33375`, `lambda = 2/3`, `S = 0.23763`, `lambda_u = 4`, and diffusivity `1.8703212744369169` for `G = kappa = 1`, agreeing with the supplied expression to roundoff. The volume ratio is exactly `4/3`. The constitutive reduction in `submission/ANSWER.md:31–37` is correct for the stated two-index equations.

## 3. Method

**Laplace–Fourier spectral Galerkin asymptotics; valid alternative.** The solver restricts displacement and hydraulic bases to the crack, eliminates the hydraulic block after that restriction, and extracts high-frequency volume and drained-tip limits. The code implements the stated Schur system, averaged Fourier-tail corrections, and finite-basis extrapolations (`submission/solve.py:13–56,101–106`). The factor `2/sqrt(pi)` converting the Laplace volume correction to the time-domain square-root coefficient is correct.

Elastic endpoint formulas and the loading-2 back-stress argument handle the simpler quantities. Resolving the tip layer is justified: the key's old face-only diffusion/reciprocity route omits part of the leading volume correction. The solver does not numerically march through the entire transient or invert a full time history. Its numerical convergence estimates are appropriately identified as nonrigorous.

## 4. Errors

**No established errors: decisive 0, side 0, presentation 0.** No incorrect reasoning was found behind the accepted numbers. The failed zero-coupling attempt was repaired in the delivered code (`submission/solve.py:48`; disclosure at `submission/ANSWER.md:123`). Unresolved high-frequency tip fits and inadequate-cutoff offsets are disclosed and excluded from the final intensity estimate (`submission/ANSWER.md:99,231`); they are informative failed checks, not concealed final-answer errors. The checkable challenges are not counted as errors. Failure to name a source is recorded separately, since source attribution was not a solver deliverable in `TASK_PROBLEM.md`.

## 5. Checks

Count **16 reported comparisons**, separating comparisons within shared headings and counting the expressly reported final reproduction comparison again. A paired volume/intensity limit comparison counts as one test; individual parameter settings within a convergence comparison do not inflate the count.

| # | Reported comparison; answer lines | Quantity tested and ability to detect an error |
|---|---|---|
| 1 | Reduced versus supplied diffusivity; 113–120 | Constitutive/time scaling; can fail for an incorrect reduction. |
| 2 | Reduced versus elastic undrained Lamé constant; 113–120 | Undrained constitutive state; can fail for incorrect storage or trace coefficients. |
| 3 | Uncoupled spectral volume/intensity versus elastic limits; 113–123 | Volume and intensity normalization; can fail, although it does not exercise coupled tip terms. |
| 4 | Small-`q` spectral response versus drained endpoints; 125–138 | Loading-1 final volume and intensity; can fail for wrong symbols or normalization. |
| 5 | Basis-size comparison; 142–166 | Short-time volume coefficient; can reveal inadequate spatial resolution. |
| 6 | Cutoff, spacing, and overlapping volume extrapolations; 168–189 | Short-time volume coefficient; can reveal quadrature, tail, or extrapolation errors. |
| 7 | Initial intensity across basis sizes, cutoffs, and fits; 193–231 | The deciding actual-tip intensity; can fail and explicitly exhibits unresolved fits. |
| 8 | Production cutoff doubling; 235–254 | Both deciding short-time quantities; can fail for inadequate Fourier tails. |
| 9 | Recomputed stdout versus saved JSON; 258–262 | Reproducibility only; detects stale or inconsistent artifacts, not a consistently wrong physical answer. |
| 10 | Required keys present; 262 | Completeness only; does not test numerical correctness. |
| 11 | Integrated drained loading-1 opening versus volume formula; 263 | Final volume, factor of two, and units; can fail, though it shares the assumed elastic opening. |
| 12 | Integrated loading-2 correction versus volume formula; 264 | Closure volume and opening normalization; can fail, but does not independently establish the back-stress model. |
| 13 | Final stdout versus saved JSON; 307 | Repeated reproduction check; artifact consistency, not mechanics. |
| 14 | Final JSON block versus answer text; 308 | Artifact consistency only. |
| 15 | Deliverables exist; 309 | Artifact completeness only. |
| 16 | Every formula is a string; 310 | Required representation only. |

**10 checks could fail on a substantive error in the model, state, or computed quantities.** The other six are falsifiable artifact checks, but a consistent, complete, physically wrong solution passes them; they are excluded from `checks_could_fail`. Neither convergence checks nor opening-integral consistency establish correctness against every possible error.

## 6. Claims

**Supported: yes; no mismatches found.** `AUTO_CHECKS.json` reports successful reproduction, all 12 tracked claims printed, and all corrected-key and principle checks passing. Its pre-P374 failures concern exactly the two subsequently corrected key values.

I reran `python3 solve.py` in `scratch/reproduction/`. It exited successfully with empty stderr; stdout and the newly written JSON both exactly equal the delivered JSON, including formulas and uncertainty metadata. The answer's pasted JSON also matches. The 13 top-level JSON fields include the reproduction-command field; this is consistent with the automatic layer's 12 tracked claims.

Additional scratch runs reproduced **all 14 printed spectral tables and all 10 groups of quoted fitted values**, including the disclosed bad fits and production cutoff doubling. Production fits are `0.8075897671067087` and `2.4400664995666825`; doubled-cutoff fits are `0.807587689880131` and `2.440069858483082`, supporting the reported rounding. Independent volume integration and constitutive arithmetic also pass. Evidence is saved in `scratch/reproduction_audit.json`, `scratch/key_comparison.json`, `scratch/report_checks.json`, and `scratch/report_claims_audit.json`. Original submission file hashes are unchanged. Check drivers are not bundled as a separate test suite in the submission, but the delivered spectral function and stated configurations reproduce their results.

## 7. Questions and missing information

**Eight listed interpretation/uncertainty items:** five bullets under “Questions and Assumptions” (`submission/ANSWER.md:49–53`) and three under “Missing or Uncertain” (`submission/ANSWER.md:103–105`). There are no direct unanswered interrogative questions. The items sensibly address shear/symmetry and contact, outer versus boundary pressure, the two-dimensional long-time limit, the incompressible endpoint, the tip/time order, absent elementary formulas or rigorous error bounds in this derivation, tip-resolution sensitivity, and the distant-field assumption for loading 2. They do not prevent completion.

## 8. Source

**None.** The submission names no source paper or author for this problem. NumPy and SciPy are implementation dependencies, not source attribution. The Detournay–Cheng attribution in `KEY.md` cannot be credited to the solver.

## 9. Challenges and key assessment

1. **Unspecified initial-tip limit.** Quote: “The initial-tip limit is not specified explicitly” (`submission/ANSWER.md:53`). The argument distinguishes an outer undrained elastic field from the drained zone at the actual tip, giving different orders of spatial and temporal limits. It supplies the nested limit, a spectral computation of approximately 0.808, and the outer elastic value 1. **Checkable: yes**, by repeating the spectral high-frequency limit and comparing the outer undrained crack solution. Its warning about equating elastic energy-release rates without accounting for Darcy dissipation belongs to this argument. Against the corrected key: (a)/(c) right. Against the actual-tip challenge: (a)/(c) right; under the reverse-order reading, replace only `K1_0` by 1. These answers are marked contested.

2. **Incomplete face-only short-time approximation.** Quote: “A one-dimensional face-only diffusion approximation misses contributions of the same order as the requested volume growth” (`submission/ANSWER.md:57`). The crack-restricted spectral system retains hydraulic tip structure and converges to `A = 2.44`, rather than the old face-only `1.354055`. **Checkable: yes**, using the specified basis/cutoff/spacing comparisons and high-frequency volume extrapolation. The difference is approximately 1.086, agreeing with the key's independently recorded tip contribution. Against the corrected key and against this challenge, (a)/(c) are right; against the superseded paper coefficient they would be partly right. This challenge is not an error.

The additional cautions about boundary pressure, absent unilateral contact, and a nonexistent decaying two-dimensional steady pressure field are interpretation notes consistent with the stated transient model, rather than allegations that the current key's loading-2 answer is wrong. The steady-field caution is checkable: a nonzero net flux in a two-dimensional harmonic exterior produces a logarithmic far field, incompatible with decay to zero; the relevant limit is instead taken from the transient. These notes do not create a conflicting answer to (b).

**Key assessment:** the active corrected numerical targets are consistent with the submission and checks. The old `0.75` initial intensity and `1.354055` coefficient retained in the historical section must not be used as current grading targets. The submission supports the key's documented correction.

There is also a dimensional typo in the historical opening–intensity relation at `KEY.md:61`: its prefactor is printed as `sqrt(pi*G)/(sqrt(8)*(1-nu))`, whereas it must be `G*sqrt(pi)/(sqrt(8)*(1-nu))`. Substituting the key's drained ellipse gives `lim(w/sqrt(xi)) = 2*sqrt(2)*(1-nu)*sqrt(L)/G`; the corrected prefactor returns the stated `K = sqrt(pi*L)`, while the printed prefactor spuriously divides it by `sqrt(G)`. This is a key presentation error, not a solver error. At the unit-modulus normalization used for the numerical comparison the two prefactors coincide; all answer grades are unchanged under either numerical reading.

```json
{"answers": {"a": "contested", "b": "right", "c": "contested"},
 "supporting_state": "right",
 "method": "Laplace-Fourier spectral Galerkin asymptotics", "method_class": "valid_alternative",
 "errors": {"decisive": 0, "side": 0, "presentation": 0},
 "checks_reported": 16, "checks_could_fail": 10,
 "claims_supported": "yes",
 "questions_listed": 8,
 "source_named": "none",
 "challenges": [{"claim": "Initial-tip intensity depends on the unspecified order of spatial and time limits", "checkable": true}, {"claim": "Face-only diffusion omits a leading tip contribution to short-time volume growth", "checkable": true}]}
```
