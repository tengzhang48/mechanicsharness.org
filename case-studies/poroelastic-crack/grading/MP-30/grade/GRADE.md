# Submission grade

## 1. Answers

Use the corrected P374 key, with 1% tolerance for nonzero values and absolute tolerance 0.01 for `K2_0`. The superseded paper values 0.75 and 1.354055 are not the grading targets.

| Part | Grade | Reason |
|---|---|---|
| (a) | **Partly right** | `K1_inf` is right; `K1_0 = 1` receives the key's outer-limit partial credit, but both volumes are doubled and the small-time coefficient is wrong. |
| (b) | **Partly right** | `K2_0 = 0` is right numerically, although the stated traction-only argument is invalid; both long-time answers are wrong, including the sign of the volume. |
| (c) | **Partly right** | Two of nine numerical values pass; the initial loading-1 intensity additionally has the limited interpretation allowed by the key. |

The numbers are in the requested units:

| Quantity | Submitted | Corrected key | Assessment |
|---|---:|---:|---|
| `K1_0` | 1 | 0.807593 | Partly: correct undrained outer intensity, not the tip-first intensity |
| `K1_inf` | 1 | 1 | Right |
| `V1_0` | 3.7699111843 | 1.8849555922 | Wrong: factor 2 |
| `V1_inf` | 5.0265482457 | 2.5132741229 | Wrong: factor 2 |
| `V1_short_coeff` | 1.4179630807 | 2.440220 | Wrong: outside the band; not the complete short-time result |
| `p_faces_0` | 1.2484394507 | 0.7490636704 | Wrong: extra factor `1/(1−ν_u)` |
| `K2_0` | 0 | 0 | Right number, defective argument |
| `K2_inf` | 0 | −0.33375 | Wrong |
| `V2_inf` | 4.4736279387 | −0.8388052385 | Wrong magnitude and sign |

For `K1_0`, reversing the limit order makes the submitted 1 correct as an outer elastic intensity. With the key's order, taking the crack-tip limit at each positive time before taking time to zero, it is incorrect. Parts (a) and (c) remain partly right under either reading because their other errors persist. The submission does not distinguish these limits.

## 2. Supporting state

**Wrong for the requested supporting quantity.** The requested `p_faces_0` is incorrect. In the stated two-dimensional equations, `ζ = 0` gives

`p = −(ν_u−ν) σ_kk/[α(1−2ν)]`.

With the correctly identified outer face stresses `σ_xx = σ_yy = −1`, this gives 0.7490636704. Inserting the submitted pressure instead leaves the right side of the given `2Gζ` equation equal to **0.712**, rather than zero. The calculation in `derive.py:120` introduces `σ_zz` into a trace explicitly restricted to indices 1 and 2.

Credit the correct two-dimensional effective-stress identity, undrained constitutive reduction, outer in-plane face stresses, and local half-space diffusion profile. The supplied opening is nevertheless doubled, and the extra out-of-plane calculation does not describe the stated reduced equations.

## 3. Method

**Elastic effective-stress reduction with ad hoc one-dimensional pressure averaging — general_or_wrong.** The undrained reduction and the idea of face diffusion are appropriate starting points. The solver then treats the total-stress intensity as fixed solely by face traction and treats effective stress as an equilibrated ordinary elastic field. In fact, total-stress equilibrium implies `σ^eff_ij,j = α p_,i`, so the effective-stress transformation alone cannot solve the pressure-loaded elastic boundary problem.

The volume calculation also starts from an incorrect elastic compliance. Its small-time calculation assumes a flat pressure influence over depth `L`, rather than deriving the mechanical weight through reciprocity and accounting for the tip contribution required by the corrected key. Neither a valid coupled solution nor a valid alternative asymptotic derivation is delivered. The correct `K2_0` number is not justified by the claim that zero face traction prevents a stress intensity at all times.

## 4. Errors

Count distinct issues below once, grouping their repeated appearances in the answer, formulas, and code. Explicitly tentative scratch expressions are not treated as additional asserted results.

### Decisive: 6

1. **Wrong instantaneous tip limit.** `ANSWER.md:68–70`: “`K1_0 = K1_inf = 1`” and “the undrained/drained distinction does not change K.” `solve.py:53` implements this. It omits the drained tip zone and confuses the outer undrained intensity with the requested tip-first limit; the partial-credit exception above applies.
2. **Doubled elastic opening and volume.** `ANSWER.md:76–85`, `check_elastic.py:27`, and `solve.py:54–55` use `w = 4(1−ν)P/G sqrt(L²−x²)`. The total opening is `2(1−ν)P/G sqrt(L²−x²)`, and its volume is `π(1−ν)PL²/G`. Independently, the plane-strain traction integral for `w=B sqrt(L²−x²)` gives `σ_yy(0)=−GB/[2(1−ν)]`; their prefactor produces **−2P**, not the prescribed −P. This changes both loading-1 volumes and also contaminates their loading-2 formula.
3. **Inconsistent three-dimensional pressure trace.** `derive.py:120–124` sets `σ_kk = −2 + σ_zz`, despite the stipulated two-dimensional trace. Consequently `ANSWER.md:116` and `solve.py:57` include an erroneous `1/(1−ν_u)`. This changes the requested initial pressure and the claimed supporting out-of-plane state.
4. **Unjustified small-time mechanical influence.** `ANSWER.md:97–101` assumes `S=(1/L)∫_0^L p/p0 dy` and infers `A=(V1_inf−V1_0)2/sqrt(π)`. Choosing the only macroscopic length does not determine an influence kernel. The resulting 1.417963 is neither the corrected 2.440220 nor a derivation of the missing tip contribution. Fixing only the doubled volumes would give 0.708982, still wrong.
5. **Wrong loading-2 long-time intensity.** `ANSWER.md:71–72,129` and `solve.py:61` assert zero total face traction means `K2_inf=0`. Pressure-induced stress throughout the solid can produce a tip intensity with zero face traction. The correct long-time intensity is `−η sqrt(πL)`, where `η=α(1−2ν)/[2(1−ν)]`.
6. **Wrong loading-2 volume mechanism and sign.** `ANSWER.md:123–127` says the crack “opens under the effective face tension” `+α`; `verify_effective_stress.py:19–26` asserts an ordinary elastic crack reduction. It omits the effective-stress body-force term and poroelastic back stress, and treats positive face normal stress as an opening pressure despite loading 1's negative normal-stress convention. The long-time result is a closure, `V2_inf=−απ(1−2ν)L²/(2G)`, not `+α V1_inf`.

### Side: 9

1. **False independent confirmation of the opening prefactor.** `ANSWER.md:133–143,176–178` claims the quadratures could expose choosing 2 instead of 4. In `check_elastic.py:27,31,81–83` and `verify_V1_0.py:40–42`, both sides share the assumed prefactor; changing it consistently preserves agreement. The purported independent opening solver raises `NotImplementedError` at `check_elastic.py:60`. These checks establish integration consistency, not the asserted elastic solution.
2. **False claim of a coupled/exact poroelastic validation.** `diffusion.py:3–4` says “verified against the exact poroelastic crack solution”; `check_A.py:5,14–15` says it solves the “coupled” problem “WITHOUT assuming the influence weight.” The code only integrates the assumed flat average over `[0,L]`. There is no independent coupled mechanical solution or weight calculation.
3. **Unsupported numerical convergence claim.** `ANSWER.md:103–104` says the numerical slope “converged” as time decreases to `10⁻⁸`. The fixed-grid errors actually grow from about `1.48×10⁻⁸` at `10⁻⁶` to `1.48×10⁻⁶` at `10⁻⁸`. The values are close to the assumed formula, but this is not a demonstrated convergence under refinement.
4. **Delivered diffusion check is unreachable.** `diffusion.py:113–123` places its comparison/printing block after `return erf(z)` at line 112, inside `erf_helper`. Running the delivered file exits successfully with no output, so it does not perform the reported FD comparison.
5. **FD comparison uses different coordinates.** `diffusion.py:120–123` selects the FD grid point nearest `yv` but evaluates the exact solution at `yv`, then labels both with `y[j]`. For the pasted `t=0.5` example, FD at `y=0.5941677521` is 0.5525443866; the exact deficit there is 0.5523999191, not the pasted 0.5485062355, which is evaluated at 0.6. The claimed difference includes a coordinate mismatch.
6. **Pressure derivation's numerical reproduction fails.** `derive.py:119,135–137` substitutes symbols without the positive assumptions attached to the symbols in its returned expression. The substitution does not evaluate it, and `float(val)` raises `TypeError: Cannot convert expression to float`. Also, its separate undrained-reduction `main()` is never called by the entry point. The valid zero residual is supplied by `verify_V1_0.py`, not by running this derivation file.
7. **Wrong square-root branch in the purported exact global stress field.** `mode1_stress.py:29` uses the principal root without choosing the branch asymptotic to `z`. On the negative real side its stresses approach −2P rather than zero at infinity: at `x=−100L`, `y=10⁻⁶L`, `σ_yy≈−2.000050P`; at `+100L`, it is approximately `0.000050P`. The reported checks sample only nonnegative x, so this does not alter the submitted scalar answers.
8. **Wrong unrequested instantaneous loading-2 opening.** `loading2.py:29–31` suggests an immediate undrained opening `2π(1−ν_u)αL²/G`. A pore-pressure-only step has no finite instantaneous face displacement; its immediate volume is zero. This extra assertion is not used in a requested volume answer.
9. **Wrong general conversion of a time slope.** `check_A.py:52` says “slope*sqrt(c)/L should equal A/L.” Its numerical slope is taken against `sqrt(t)` (`check_A.py:39`); the correct relation is `A = slope*L/sqrt(c)`. Setting `L=c=1` hides the error in this run.

### Presentation: 3

1. **Pressure/deficit mislabeled.** `ANSWER.md:159–164` calls the FD comparison a test of `p=p0 erf(...)`, but its boundary value 1 and interior values are the deficit `u=p0 erfc(...)`. Actual face pressure is zero. The same mislabel appears in `diffusion.py:115`.
2. **Overstated printed coverage.** `ANSWER.md:3–4` says `solve.py` “prints every value in results.json.” It prints nine scalar entries, but does not print the three additional top-level entries `formulas`, `reproduce_command`, and `parameters`. It does write them correctly.
3. **Contradictory boundary wording.** `derive.py:60–62` calls the boundary “traction-free” while assigning its normal traction −1. Loading 1 is a loaded face, with zero shear traction.

## 5. Checks

Count the **seven numbered check groups** in `ANSWER.md:131–182`, without recounting the same undrained algebra under group 6. “Could fail” below means an actual delivered, non-circular comparison capable of detecting an error in the physical/algebraic quantity it purports to check; file reproduction is assessed separately.

| Reported check | Does it test a deciding quantity? | Count as could fail? |
|---|---|---|
| 1. Elastic opening/volume quadrature and tip coefficient | No independent traction/compliance test: both sides use the same wrong opening | No |
| 2. Effective-stress and undrained symbolic reductions | Yes for the constitutive identities, though not for equilibrated crack stress or volume | **Yes**: a wrong coefficient or trace leaves a nonzero residual |
| 3. Face stress state | Yes for the outer in-plane trace entering the pressure calculation, not for the later erroneous 3-D trace | **Yes**: the evaluated stress field could disagree with the asserted face state |
| 4. FD diffusion profile | At most the local diffusion profile; it cannot determine the mechanical influence or complete small-time coefficient | No as delivered: the comparison is unreachable; its helper solves the chosen heat equation, and the pasted comparison also mismatches coordinates |
| 5. Small-time coefficient | It tests quadrature of the assumed average, not the actual poroelastic volume | No |
| 6. `V1_0` quadrature | Same assumed doubled opening; its valid algebra was already counted under group 2 | No |
| 7. Whole-solution reproducibility | File consistency, not physical correctness | No: all seven out-of-band values reproduce |

Thus **7 reported, 2 counted as capable of failing**. The FD helper is a real numerical diffusion solver and was exercised in this review; that does not repair the delivered verification driver or validate its crack-volume coupling. The two automatic principle checks pass, but were not additional solver-reported checks: a common volume factor cancels in the ratio, and `K1_inf=1` is correct.

## 6. Claims and delivered evidence

**Partly supported overall.** `AUTO_CHECKS.json` records reproduction exit 0, nine printed scalar claims out of twelve top-level entries, two passing value checks, and both passing principle checks. The independent rerun agrees: `python3 solve.py` exits 0, prints all nine submitted scalar numbers, and regenerates an object exactly equal to the delivered `results.json`, including formulas and metadata. The reproduction command does not run the seven checks.

The pasted elastic and undrained-volume quadratures, symbolic residuals, and `check_A.py` outputs reproduce. The cited face-state values agree with evaluation of its helper, although the complete pasted line is not the default script's output at that point. `derive.py` exits 1 before printing its numerical pressure; `diffusion.py` prints nothing. Their reported verification is therefore not reproduced by those entry points. The FD helper can reproduce the deficit numbers with the sampling problem described above. None of the reproducible quadratures independently verifies the wrong opening or short-time model.

Detailed entry-point stdout, stderr, exit codes, and recorded child PIDs are in `scratch/rerun_logs.json`; numerical band comparisons are in `scratch/key_comparison.json`. The original submission was left untouched.

## 7. Questions and missing information

There are **five listed topics** under `ANSWER.md:184–206`, counted as the solver's questions/assumptions list:

1. Opening versus half-opening: a convention check is reasonable generally, but the defined crack volume uses total opening, and the factor-of-two error is in their elastic formula, not missing information.
2. Influence depth/weight: a sensible admission of a genuine derivational gap. It gives no basis for selecting their coefficient; the missing information is the mechanical calculation, not additional problem data.
3. Initial face-side pressure: the order of the undrained outer limit versus the imposed drained boundary merits care. Their interpretation does not justify adding a third component to the stipulated trace.
4. Total versus effective loading-2 intensity: not a missing definition. A total-stress intensity can be generated by poroelastic coupling with zero applied face traction; their hypothetical `α sqrt(πL)` is not established.
5. Positive small-time growth: a correct qualitative sanity observation, rather than a question or missing input.

## 8. Source

`check_elastic.py:10` names **“Hills, Kelly, Dai, Korsunsky, \"Solution of Crack Problems\"”**. This is an elastic-crack reference, not identification of the poroelastic source problem. The submission also invokes Westergaard in `mode1_stress.py:7`. It does not name Detournay and Cheng or their 1991 paper. Credit the explicit elastic reference only; it does not substantiate the reported doubled prefactor.

## 9. Challenges and key observations

**No substantive solver challenge, and no contested answer.** The convention alternatives and influence-weight uncertainty in section 4 of the answer are assumptions, not a demonstrated contradiction in the problem, the key's model, or a published result. No coupled computation, limit calculation, or independently derived counterexample is offered to defend the wrong values. In particular, the partial credit for the outer `K1_0=1` comes from the key's explicit order-of-limits allowance, not a solver-argued challenge.

There is a non-decisive dimensional typo in the key's transcription of the opening–intensity relation: its Eq. 27 prefactor is written `sqrt(πG)/[sqrt(8)(1−ν)]`. With `G` a stress and `D_n/sqrt(ξ)` a square root of length, the required prefactor is `G sqrt(π)/[sqrt(8)(1−ν)]`. Correcting that transcription leaves all stated numerical bands and all answer grades above unchanged; accepting the transcription as printed also gives the same normalized comparisons at `G=1`. This is a grader observation, not a solver challenge or an error charged to the submission. No evidence in this submission warrants replacing the corrected P374 values.

```json
{"answers": {"a": "partly", "b": "partly", "c": "partly"},
 "supporting_state": "wrong",
 "method": "elastic effective-stress reduction with ad hoc 1-D averaging", "method_class": "general_or_wrong",
 "errors": {"decisive": 6, "side": 9, "presentation": 3},
 "checks_reported": 7, "checks_could_fail": 2,
 "claims_supported": "partly",
 "questions_listed": 5,
 "source_named": "Hills, Kelly, Dai, Korsunsky, \"Solution of Crack Problems\"",
 "challenges": []}
```
