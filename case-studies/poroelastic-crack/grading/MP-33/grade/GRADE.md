# Grading profile

## 1. Answers

Use the **corrected** key in `KEY.md:11–33`, rather than the paper's superseded values. Nonzero numerical answers have a 1% band; `K2_0` has absolute tolerance 0.01. K values below are normalized by √(πL), and V values by L²/G.

| Question | Judgment | Reason |
|---|---|---|
| (a) | **Partly right** | `K1_inf`, both endpoint volumes, and the induced face pressure are right; `K1_0=1` is the undrained outer limit, deserving the key's explicit partial credit, while the tip-first limit is 0.807593; the growth coefficient 1.128379 is wrong against 2.440220. |
| (b) | **Partly right** | `K2_0=0` is right, but `K2_inf=0` and `V2_inf=0` are wrong against −0.33375 and −0.838805; the proposed long-time argument is assessed in section 9. |
| (c) | **Partly right** | Five of nine reported numerical values meet the bands; four do not, although `K1_0` is correct under the reverse order of limits acknowledged by the key. |

| Quantity | Submitted | Corrected key | Assessment |
|---|---:|---:|---|
| `K1_0` | 1 | 0.807593 | Outside band by 23.82%; partial credit for outer-limit interpretation. |
| `K1_inf` | 1 | 1 | Right. |
| `V1_0` | 1.8849555922 | 1.8849555922 | Right. |
| `V1_inf` | 2.5132741229 | 2.5132741229 | Right. |
| `V1_short_coeff` | 1.1283791671 | 2.440220 | Wrong; 53.76% too small. |
| `p_faces_0` | 0.7490636704 | 0.7490636704 | Right. |
| `K2_0` | 0 | 0 | Right. |
| `K2_inf` | 0 | −0.33375 | Wrong. |
| `V2_inf` | 0 | −0.8388052385 | Wrong. |

The general formulas for the five matching values are also correct. The other four formula strings implement the same incorrect or incomplete conclusions as their numbers.

## 2. Supporting state

**Right for the intermediate state explicitly requested:** the induced, pre-diffusion face-adjacent pressure is correctly obtained as `2(ν_u−ν)/[α(1−2ν)]`, giving 0.7490636704. This is the undrained outer pressure before the prescribed zero-pressure boundary condition drains the adjacent solid. The undrained constitutive reduction is also right. Some additional, unrequested intermediate reasoning is faulty, as itemized below.

## 3. Method

**Undrained/drained elasticity plus a local 1D face-diffusion layer; class: general_or_wrong.** Elasticity correctly supplies the endpoint volumes, loading-1 long-time K, and the outer induced pressure. The needed asymptotic method must also resolve the drained crack-tip limit and use a global volume relation, including the leading tip contribution. Multiplying local face heave by crack length does not supply those pieces. For loading 2 the solver substitutes a uniform free-swelling state without establishing that it is the limit of the prescribed infinite-domain evolution. This is an incomplete/wrong route for the full requested solution, rather than a validated alternative.

## 4. Errors

Counts are **4 decisive, 8 side, 1 presentation**. Repeated appearances of an erroneous result in the answer, code, and JSON are counted once; separately false assertions are listed separately.

**Decisive**

1. **Neglecting the tip limit.** `submission/ANSWER.md:17–18`: a zero-thickness pressure layer “cannot affect K or V at t = 0⁺”; lines 28–30 consequently set `K1_0=1`. Its vanishing size does not make its effect on a singular tip intensity vanish. The key's tip-first value is 0.807593. The reported 1 remains correct for the reverse-order outer limit, so this receives partial credit rather than being wholly dismissed.
2. **Local heave is substituted for the global crack-volume response.** `submission/ANSWER.md:48–50`: “Both faces over length 2L give ΔV = 4LΔu.” This omits the elastic adjustment captured by reciprocity. Even the key's face-only contribution is `16(1−ν_u)(ν_u−ν)/[(1−ν)√π] = 1.354055`, whereas this construction gives 1.128379, a factor 1.2 smaller for these parameters.
3. **Tip contribution assigned the wrong asymptotic order.** `submission/ANSWER.md:50–51`: “tip zones of size √(ct) contribute only O(ct/L²).” The corrected key assigns a further **1.086165 to the coefficient of √(ct/L²)**. Omitting this leading contribution persists even if the face-only calculation is repaired. The same false order assertion recurs at lines 120–121.
4. **Pointwise pressure convergence is mistaken for zero limiting total stress.** `submission/ANSWER.md:61–65`: “Then σ → 0 … uniform free swelling … leaves the crack invisible.” The expanding pressurized region retains a local confining stress −η, with `η=α(1−2ν)/[2(1−ν)]`. Canceling its normal traction on the crack gives `K2_inf=−η` and `V2_inf=−απ(1−2ν)/2` in the reporting units, rather than zero. See the checkable counterexample in section 9.

**Side**

1. **Total and incremental face-layer strains are conflated.** `submission/ANSWER.md:43–47` writes `e_xx=0` and total-stress/strain formulas with no baseline term. The undrained face state actually has `e_xx=−(1−2ν_u)P/(2G)=−0.1P/G`. Substitution of its `p0` into the written trace formula gives −1.75P instead of −2P. The relations work for increments with `δe_xx=0`; the missing constant cancels from the pressure derivative and heave difference, so it does not change the reported diffusivity or numerical coefficient.
2. **Wrong square-root branch in the purported Westergaard field.** `submission/solve.py:62` uses the principal `np.sqrt(z*z−L*L)`, rather than the branch asymptotic to z. On the negative x axis its Z tends to −2P, violating decay at infinity; at x=−100L the real value is −2.000050P. Its limiting interior-face real trace still gives −2P, so the reported face pressure survives this error.
3. **A tautology is described as a diffusivity consistency test.** `submission/solve.py:44–47` compares two algebraically identical copies of the stated c; `submission/ANSWER.md:90–91` says it fails if the mechanical constraint is wrong. No derivation from storage or a changed constraint is evaluated by this assertion.
4. **The face-pressure check does not use the computed face trace.** `submission/ANSWER.md:99–100` says pressure is “recomputed from that trace.” In `submission/solve.py:67–68`, the trace is hard-coded as −2P, reproducing the same formula as `p0`; `skk_face` is unused in this pressure assertion.
5. **Overstated physical coverage of the numerical layer checks.** `submission/ANSWER.md:107–108` says the FDM test fails if “any element of the layer physics (diffusivity, heave law, both-faces factor 4L, erfc profile)” is wrong. The simulation sets c=1 and uses the same assumed heave law and factor 4L as the analytic coefficient. It can detect numerical/algebraic disagreement within that model, but cannot validate those common mechanical assumptions or detect the missing tip response. Likewise, the COD and erfc quadratures validate integrals, not the physical validity of their supplied fields.
6. **Incorrect explanation of the first-order convergence.** `submission/ANSWER.md:105–106` and `submission/solve.py:99` attribute it to the discontinuous initial condition. The right-endpoint sum at `solve.py:86` omits the boundary half-cell term `p0*dy/2`. Adding that term changes error ratios from approximately 0.50 to 0.25 on refinement. The first-order Richardson cancellation nevertheless gives the correct coefficient of the solver's chosen model.
7. **False combined-loading endpoint claim.** `submission/ANSWER.md:66–67` asserts `K1_0+K2_0=K1_inf+K2_inf=1` for a real fluid-pressurized crack. The corrected tip-first endpoints are approximately 0.807593 and 0.66625, respectively. Superposition is valid; the asserted endpoint values are not.
8. **Incorrect description of a finite pressure sink.** `submission/ANSWER.md:132–134` says such a sink “contradicts ‘vanish at infinity’.” A finite sink changes the specified boundary-value problem; a zero-pressure sink does not, by itself, contradict vanishing pressure at infinity. This does not rescue the submitted infinite-domain limit.

**Presentation**

1. `submission/solve.py:74` says “run three grids”; the code actually runs five at lines 92–93. The answer and reproduced output correctly report five.

## 5. Checks

The answer groups its tests into seven numbered items, but the delivered code contains **11 pass/fail assertions**. Auxiliary slope/ratio printouts are diagnostics, not additional assertions. All eleven pass on rerun.

| Assertion | What it actually tests; deciding-quantity coverage | Could fail? |
|---|---|---|
| Undrained reduction | Constitutive algebra at ζ=0; does not test the tip intensity. | Yes, on an algebra/implementation error. |
| 1D diffusivity | Identical formulas for c; no independently obtained storage coefficient. | No. |
| COD volume quadrature | Integral of the supplied drained COD, relevant to `V1_inf`; common prefactor errors would survive. | Yes. |
| Integral of erfc | Elementary integral used in the layer calculation; does not test full crack-volume mechanics. | Yes. |
| Westergaard face trace | Interior-face elastic trace, relevant to `p_faces_0`; does not test the tip limit or global branch. | Yes. |
| Pressure “from trace” | Same formula as `p0`, with hard-coded trace; no independent pressure test. | No. |
| Finest-grid FDM slope | Numerical slope versus analytic A **within the assumed 1D heave model**; not the full deciding growth coefficient. | Yes, for disagreement within that model. |
| Richardson FDM slope | Convergence to that same model's analytic coefficient; not the complete crack asymptotics. | Yes. |
| erf-profile integral, ct=10⁻⁴ | Numerical integral versus its analytic integral; local profile only. | Yes. |
| erf-profile integral, ct=10⁻³ | Same test at another scale. | Yes. |
| erf-profile integral, ct=10⁻² | Same test at another scale. | Yes. |

**Checks that could fail: 9.** This counts non-tautological assertions on their stated numerical/algebraic targets; it does not mean nine independent checks of the complete physical solution. No reported check evaluates the drained-tip K, the loading-2 limiting confinement, or the full leading volume coefficient. The two principle relations checked by `AUTO_CHECKS.json` also pass: the endpoint volume ratio and the Griffith long-time `K1_inf`.

## 6. Claims and reproduction

**Numerical/output claims supported: yes.** `AUTO_CHECKS.json` records successful reproduction and nine of nine scalar claims printed. A fresh `python3 solve.py` run on the copy in `scratch/reproduction/` exits successfully, prints all nine stated values, reports exactly `11 checks, 0 failed`, and generates a `results.json` byte-for-byte identical to the submission. The pasted check values, grid slopes, ratios, and extrapolation agree with the rerun. There are no numerical reproduction mismatches. Claims about what the checks establish are overstated as documented in sections 4–5; reproducing a number does not establish its correctness.

## 7. Questions and missing information

**Five listed assumption/question items**, at `submission/ANSWER.md:125–137`: loading magnitude/units, sign convention, loading-2 long-time pressure conditions, quasi-static response, and face permeability. They contain no unresolved request for information. Four restate conventions or conditions already provided. The two-dimensional far-field limit is a sensible issue to investigate, but its proposed mechanical conclusion is wrong. The separate three uncertainty bullets concern unrequested intermediate-time behavior and higher-order corrections; they are not additional questions or missing data needed for this task. Their assertion that tip corrections are higher order is already counted above.

## 8. Source

**None identified for the poroelastic crack problem.** The answer mentions the “Westergaard solution” (`submission/ANSWER.md:28`) and “Terzaghi effective stress in 1-D” (line 47), along with Griffith crack terminology. These name tools or eponymous models, not the source paper or an author attribution for this problem. It does not name Detournay and Cheng or their 1991 paper.

## 9. Challenges and key assessment

**One checkable implicit challenge:** loading 2 cannot have a decaying two-dimensional steady pressure field, so its pressure becomes locally uniform and its long-time crack response must vanish (`submission/ANSWER.md:59–67,130–134`). The first observation is sound: a bounded harmonic pressure with constant value 1 on the finite crack cannot also decay at infinity. The proposed jump from pointwise p→1 to σ→0 is not justified, and no delivered calculation tests it: the code simply assigns the zero endpoints.

A mechanical counterexample is reproducible in `scratch/audit.py`. Take a circular pressure patch p=1 for r<R and p=0 outside, in the same infinite plane-strain elastic medium. Let

`a = α(1−2ν)/[4G(1−ν)]`,

with radial displacement `u_r=a*r` inside and `u_r=a*R²/r` outside. These fields satisfy equilibrium, displacement and radial-traction continuity, and vanishing stress/pressure at infinity for every finite R. Nevertheless, the interior stress is

`σ_rr=σ_θθ=−α(1−2ν)/[2(1−ν)]=−η`,

independent of R. As R grows, pressure tends to 1 at every fixed point **without** the stress tending to zero. For the stated parameters the calculation gives −0.33375, not zero. This checks the disputed mechanical inference; it is not a simulation of the crack transient. The same retained confinement, followed by the elastic correction to make the crack faces traction-free, gives the key's negative loading-2 K and V.

**Against the key:** the two loading-2 long-time answers are wrong and (b) is partly right. **Under the solver's assumed uniform free-swelling state:** its zero crack perturbation is internally consistent, but that state has not been obtained as the specified evolution's limit and its determining inference is contradicted above. Thus this is a checkable but rejected challenge, not a well-argued alternative eligible for a contested grade or exemption from the decisive-error count.

For `K1_0`, the key already acknowledges the statement's order-of-limits ambiguity: 1 is right if time tends to zero before taking the tip limit, while 0.807593 is right for the tip-first limit. The solver does not analyze these competing limits and instead declares the boundary layer unable to affect K. Its outer interpretation receives the explicit partial credit required by the key. No additional error in the corrected key is established here; the older 0.75 and 1.354055 entries are explicitly superseded.

```json
{"answers": {"a": "partly", "b": "partly", "c": "partly"},
 "supporting_state": "right",
 "method": "Undrained/drained elasticity plus local 1D face diffusion", "method_class": "general_or_wrong",
 "errors": {"decisive": 4, "side": 8, "presentation": 1},
 "checks_reported": 11, "checks_could_fail": 9,
 "claims_supported": "yes",
 "questions_listed": 5,
 "source_named": "none",
 "challenges": [{"claim": "Loading 2 has zero long-time crack response because pressure becomes locally uniform", "checkable": true}]}
```
