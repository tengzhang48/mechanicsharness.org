# Grade

## 1. Answers

Use the corrected key in `KEY.md:11–34`, not the historical paper values reproduced later in that file. The bands are 1% for nonzero quantities and absolute 0.01 for `K2_0`.

- **(a) Partly right.** Both volumes and the long-time intensity are right; the initial intensity is the undrained outer limit rather than the drained-tip limit, and the leading volume-growth coefficient is wrong.
- **(b) Partly right.** `K2_0 = 0` is right, but both long-time quantities incorrectly vanish instead of giving the negative back-stress intensity and closure.
- **(c) Partly right.** Four of nine numbers pass; `K1_0` receives the key's stipulated partial credit for the alternative order of limits, while the other four failing quantities are wrong.

All values below use the problem's normalization.

| Quantity | Submission | Corrected key | Judgment |
|---|---:|---:|---|
| `K1_0` | 1 | 0.807593 | Partly: reverse-order, undrained outer limit only |
| `K1_inf` | 1 | 1 | Right |
| `V1_0` | 1.884955592154 | 1.884955592154 | Right |
| `V1_inf` | 2.513274122872 | 2.513274122872 | Right |
| `V1_short_coeff` | 0.564189583548 | 2.440220 | Wrong: 76.88% low |
| `p_faces_0` | 0.374531835206 | 0.749063670412 | Wrong: half the required value |
| `K2_0` | 0 | 0 | Right |
| `K2_inf` | 0 | −0.33375 | Wrong |
| `V2_inf` | 0 | −0.838805238508 | Wrong |

The general formula strings repeat these same successes and failures. The two automatic principle checks, the volume ratio and the long-time Griffith intensity, pass.

## 2. Supporting state

**Wrong for the explicitly requested auxiliary state:** the undrained pressure immediately adjacent to the faces. Eliminating fluid content correctly gives
`p = −(ν_u−ν)σ_kk/[α(1−2ν)]`, but the zero-remote-stress elastic crack has `σ_xx = σ_yy = −σ₀` on its faces, hence `σ_kk = −2σ₀`. The submitted trace is half as large.

Some additional, unrequested supporting calculations are right: the drained and undrained opening ellipses, the local incremental strain–pressure slope, and the quoted diffusivity. The claimed zero-stress long-time loading-2 background is wrong.

## 3. Method

**Undrained/drained Griffith endpoints with local 1-D oedometric drainage; class: general_or_wrong.** Elastic endpoint formulas correctly determine the volumes and loading-1 long-time intensity. They cannot determine the initial tip intensity by substituting undrained moduli everywhere. The local diffusion profile is appropriate away from the tips, but its conversion into total opening omits the reciprocal elastic response and the leading tip contribution. The 1-D simulation is a check of that local model, not a full transient solution or a valid alternative derivation of the requested coefficient.

## 4. Errors

Counts below are for distinct errors; repeated formulas and consequences in the answer, code, and JSON are grouped together.

1. **Decisive — treating the initial tip as undrained.** `submission/ANSWER.md:31`: “K1_0 = σ₀√(πL)”; implemented at `submission/solve.py:45`. At fixed positive time the tip limit must include drainage; the corrected initial value is 0.807593. The statement at `ANSWER.md:120–121` that the near-tip correction “does not affect any reported item” repeats this omission. Value 1 is nevertheless partly correct under the reverse order of limits, as the key explicitly requires.
2. **Decisive — wrong tangential face stress and stress trace.** `submission/ANSWER.md:34`: “the faces have σ_xx = 0, σ_yy = −σ₀, so σ_kk = −σ₀”; `solve.py:50` uses the resulting pressure. Both face stresses are −σ₀ for the elastic crack with stresses vanishing at infinity. This halves `p_faces_0` and also lowers the short-time coefficient. The fluid-content elimination itself is correct.
3. **Decisive — converting local layer strain directly into the entire crack-volume change.** `submission/ANSWER.md:45–48`: “Δu_y(0⁺,t) = −∫₀^∞ Δe_yy dy” followed by “ΔV = 2L·2Δu_y.” Extending the local oedometric approximation through the whole displacement integral drops the global elastic adjustment. Correcting the pressure alone gives `A = 1.128379167096` under this rule, still short of the key's face-only reciprocal contribution, `1.354055000515`. Thus the erroneous coefficient is not explained solely by the missing pressure factor.
4. **Decisive — dismissing the tip contribution to the leading volume growth.** `submission/ANSWER.md:50–51`: “The crack-tip regions … contribute only O((ct)^{3/4}L^{1/2}) to V … so the uniform face contribution … is the leading term.” An estimate of the volume confined to the inner region does not include its effect on opening elsewhere. The corrected key includes a tip contribution `1.086165` at the same `√(ct/L²)` order, raising the face contribution to `2.440220`. The boundary-only loading decomposition at `ANSWER.md:21–24` likewise does not account for relaxation of the singular initial tip pressure.
5. **Decisive — inferring zero total stress from pointwise uniform pressure.** `submission/ANSWER.md:60–63`: “Uniform p means … zero total stress” and “K2_inf = 0, V2_inf = 0”; hard-coded at `solve.py:52–53`. The expanding pressure field retains the compressive background `−η`, where `η = α(1−2ν)/[2(1−ν)]`. The traction-free crack consequently has `K2_inf = −η` and `V2_inf = −π(1−ν)η` in normalized units. The asserted return of the negative histories to zero at `ANSWER.md:64–65` is also wrong. The combined hydraulic-crack endpoint claims at `ANSWER.md:66–67` inherit errors 1 and 5: its key-based intensities are initially 0.807593 and finally 0.66625, and its final volume is 1.674468884364 rather than 2.513274122872.
6. **Side — imposing zero absolute tangential strain.** `submission/ANSWER.md:40–41`: “the field is locally y-dependent only, so e_xx = 0.” The local drainage increment can have `Δe_xx ≈ 0`; the pre-existing strain need not vanish. In fact the undrained face has `e_xx = −(1−2ν_u)σ₀/(2G)`. The local pressure derivative and diffusivity survive this error because the constant baseline strain drops out.
7. **Presentation — overstating the scope of verification.** `submission/ANSWER.md:104`: “Independent physics check of the short-time coefficient.” That calculation checks the solver's chosen 1-D model with its chosen pressure and volume conversion. It does not independently test the physical coefficient for the poroelastic crack. This misdescription does not introduce a further numerical error.

**Totals: 5 decisive, 1 side, 1 presentation.** No additional error is assigned merely for omitting a source or for unverified claims about higher-order remainders.

## 5. Checks

Count the six numbered check families in `submission/ANSWER.md:89–113`; they produce nine individual PASS comparisons. “Could fail” here means that the independently evaluated numerical side could disagree with the claimed formula for an algebraic or numerical error in the relation tested. It does not mean that the check establishes the whole physical model.

| Reported check | What it actually tests | Could fail? |
|---|---|---|
| 1. COD quadrature, drained and undrained | Correct endpoint-volume formulas from the assumed ellipses; no test of tip intensity. The code re-evaluates the formulas instead of comparing the exported volume variables. | Yes, for an ellipse-integration/formula mismatch |
| 2. `∫ erfc = 1/√π` | A mathematical integral used in the layer calculation; not the physical growth coefficient | Yes, for an incorrect integral value |
| 3. Oedometric pressure slope | The local incremental constitutive coupling; not its global applicability to crack opening | Yes, for an incorrect slope formula |
| 4. Diffusivity | The local constitutive/mass-balance reduction against the supplied diffusivity | Yes, for an inconsistent reduction |
| 5. FD diffusion, three times | Reported `A` against a 1-D diffusion model sharing the wrong `p_b` and local volume rule; excludes tips | Yes, for a numerical or amplitude mismatch within that model |
| 6. Analytic-profile quadrature | The same local volume integral against `A`; shares the preceding assumptions | Yes, for an integral/amplitude mismatch within that model |

**Six checks reported and six capable of failure in those limited senses** (nine primitive comparisons). None independently tests the wrong face pressure, the initial drained-tip intensity, the loading-2 long-time state, or the complete short-time coefficient. Checks 5 and 6 passing therefore do not resolve the decisive errors.

## 6. Claims and reproduction

**Claims supported: yes, as numerical reproduction.** `AUTO_CHECKS.json` reports exit 0 and all nine numerical claims printed. I reran `python3 solve.py` on a copy in `scratch/reproduction/`, recording its PID. The rerun also exited 0; its JSON content exactly equals the submitted JSON, and its complete text output exactly equals `submission/run_output.txt`. All nine PASS comparisons and the zero-failure summary reproduce. There are no numerical or output mismatches. Reproduction does not establish correctness against the key; the scope overstatement is recorded separately above.

## 7. Questions and missing information

**Four items listed** at `submission/ANSWER.md:132–141`, mostly interpretations rather than requests for missing data:

1. Pressure units: the chosen normalization by applied stress is sensible; the task already specifies per-unit loading magnitude.
2. Maintained loading-2 reservoir: already specified by `H(t)`. Pointwise pressure equilibration is a sensible observation, but its zero-stress conclusion is erroneous.
3. Hydraulic connection of faces: consistent with the explicit pressure boundary conditions; no missing information is identified.
4. “Immediately next to” versus the exact drained face: a sensible distinction between the undrained outer pressure and the prescribed boundary value, although the outer pressure was calculated incorrectly.

The actual ambiguity identified in the key—the order of the time and tip limits for `K1_0`—is not explicitly raised or analyzed by the solver.

## 8. Source naming

**No source paper or author is identified for this poroelastic problem.** The phrases “Griffith crack formulas” (`submission/ANSWER.md:17`) and “Westergaard function” (`ANSWER.md:33–34`) name elasticity constructions. They do not identify Detournay and Cheng (1991), or another source for the poroelastic asymptotics. Record `source_named = "none"`.

## 9. Challenges and key assessment

There is no explicit critique of a named paper. There is an **implicit, checkable challenge** to the nonzero long-time loading-2 state: `submission/ANSWER.md:59–63` argues that 2-D diffusion has no decaying steady pressure solution, so pressure tends pointwise to the reservoir value, and then assumes a freely expanding, zero-total-stress state. The pointwise pressure argument is sound; the implication about stress is not.

A concrete check is provided in `scratch/review_calcs.py` and its saved output. For a prescribed pressure `P` inside a disk of radius `R` and zero outside, take the equilibrium radial displacement `u = a r` inside and `u = b/r` outside. Displacement and traction continuity give
`a = αP(1−2ν)/[4G(1−ν)]`, `b = aR²`, and an inner total stress `−ηP`, independent of `R`. Every finite-radius state has zero stress at infinity. As `R → ∞`, pressure becomes pointwise uniform but the compression remains. This is an equilibrium counterexample to the solver's implication, not a simulation of the crack transient. Adding a traction-free crack to that compressive background gives precisely the key's negative intensity and volume.

**Against the key:** (b) is partly right. **Against the assumed free uniform-expansion state:** its zero long-time values would be right for that different mechanical state, but the solver has not shown that the stated infinite-domain loading approaches it. The challenge is checkable but does not meet the well-argued standard for a contested answer; its failed inference is counted as error 5. The rerun only solves 1-D diffusion and does not verify this challenge.

The key itself acknowledges the initial-intensity order-of-limits ambiguity, and partial credit for value 1 follows its explicit instruction. I find no additional demonstrated error in the supplied corrected key; its superseded historical formulas are not used for grading.

```json
{"answers": {"a": "partly", "b": "partly", "c": "partly"},
 "supporting_state": "wrong",
 "method": "Undrained/drained Griffith endpoints with local 1-D oedometric drainage",
 "method_class": "general_or_wrong",
 "errors": {"decisive": 5, "side": 1, "presentation": 1},
 "checks_reported": 6, "checks_could_fail": 6,
 "claims_supported": "yes",
 "questions_listed": 4,
 "source_named": "none",
 "challenges": [{"claim": "Pointwise uniform loading-2 pressure implies zero long-time intensity and volume", "checkable": true}]}
```
