# Grading profile

## 1. Answers

Use the **corrected** key at `KEY.md:11–34`, not its superseded paper values. Nonzero quantities have a 1% band; `K2_0` has an absolute 0.01 band.

- **(a) Partly right.** Both endpoint volumes and the long-time intensity are right; `K1_0 = 1` is the undrained outer limit, explicitly given partial credit by the key, but misses the drained tip limit 0.807593; the small-time coefficient is wrong.
- **(b) Contested.** Against the key: partly right (`K2_0` is right, both long-time quantities are wrong). Against the submission's proposed stress-free uniform-pressure end state: all three values are right conditionally; that end state is not established for the stated initial/far-field problem. See section 9.
- **(c) Contested.** Against the key: partly right, with five of nine values within band. Under the loading-2 challenge, its two disputed zeros are additionally accepted conditionally, but the loading-1 tip intensity and short-time coefficient remain deficient.

| Quantity | Submission | Corrected key | Assessment |
|---|---:|---:|---|
| `K1_0` | 1 | 0.807593 | Outside band; partial credit for reverse order of limits |
| `K1_inf` | 1 | 1 | Right |
| `V1_0` | 1.8849555922 | 1.8849555922 | Right |
| `V1_inf` | 2.5132741229 | 2.5132741229 | Right |
| `V1_short_coeff` | 1.8054066674 | 2.440220 | Wrong |
| `p_faces_0` | 0.7490636704 | 0.7490636704 | Right |
| `K2_0` | 0 | 0 | Right |
| `K2_inf` | 0 | −0.33375 | Contested; outside key band |
| `V2_inf` | 0 | −0.8388052385 | Contested; outside key band |

The general expressions for the five passing quantities also agree with the key. The code directly assigns the other four values from the disputed formulas; it does not independently solve for them.

## 2. Supporting state

**Right for the explicitly requested intermediate quantity.** Setting `zeta = 0` gives the correct undrained pressure/stress relation, and the reference face trace is −2 per unit traction. Thus `p_faces_0 = 2(nu_u−nu)/[alpha(1−2nu)]` is correct, with a positive sign and the intended pre-diffusion interpretation (`submission/ANSWER.md:12–14,40–43`). It is distinct from the exactly imposed face pressure zero for every positive time. The drained and undrained Sneddon volume formulas are also correct. Incorrect extensions of these states are assessed below.

## 3. Method

**Elastic endpoint limits plus uncoupled pressure diffusion and Betti reciprocity; method class: `general_or_wrong`.** Elastic asymptotics and reciprocity are appropriate tools, and the elimination yielding diffusion of fluid content is correct. The actual reduction, however, substitutes undrained elasticity for the tip limit and a scalar pressure heat equation for the coupled mechanical/diffusive response. Neither substitution supplies the needed solution. The 2-D FTCS calculation is a check of that scalar surrogate, not a transient Biot solution. The loading-2 endpoint uses the separately assessed challenge.

## 4. Errors

Counts below distinguish separate errors in reasoning or computation; repeated formulas and their numerical consequences are not counted again. The loading-2 challenge is excluded under the challenge rule.

**Decisive: 3.**

1. **Missing drained tip limit.** “`K1_0 = sqrt(pi L)`” (`submission/ANSWER.md:34–35`; `submission/solve.py:24`). This gives the outer undrained intensity, not the key's tip-first limit. The submission never distinguishes the two orders of limits. Partial credit for this number does not remove the answer-changing omission.
2. **Pressure diffusion substituted for the coupled problem.** “the tangential-strain coupling is smaller” and “`p^II = -p_faces_0 erfc(...)`” (`submission/ANSWER.md:47–58`; `submission/solve.py:59–61,83–88`). Diffusion of `zeta` does not justify identifying the entire pressure perturbation with this heat solution: pressure also contains a mechanically determined harmonic component. The resulting coefficient 1.805407 is already different from the key's face contribution 1.354055, before including tips.
3. **Leading tip contribution discarded.** “Tip zones ... contribute only `O((ct)^{3/4})` and are subleading” (`submission/ANSWER.md:59,114–115`). That estimate applies to bounded pressure in the surrogate. In the coupled crack problem, singular pressure and the singular reciprocal weight make a tip contribution of order `sqrt(ct)`. The corrected key adds 1.086165 at that order, giving 2.440220.

**Side: 6.**

1. **False whole-transient claim.** “in fact it is constant for all t for this loading” (`submission/ANSWER.md:38–39`; repeated in the `K1_0` formula). The tip intensity has different short- and long-time limits. This additional claim changes no separately requested intermediate-time answer.
2. **Displayed multiplication loses a factor two.** The product at `submission/ANSWER.md:56` evaluates to coefficient 0.902703 for `L = G = 1`; the next line asserts 1.805407. The missing factor is the magnitude two of the reference stress trace. The implemented coefficient includes that factor, so this is a derivation error without an additional numerical change.
3. **Incorrect tail for the integration domain.** `submission/solve.py:43–49` integrates over a square `[-4,4]^2` but adds the leading tail for the exterior of a circle of radius 4. Its tail 0.098175 should be approximately 0.090828 for this square and field (local quadrature in `scratch/audit.py`). This affects CHECK 1, not a reported answer. The comment's use of `|Z|` also omits the angular averaging needed for `(2 Re Z)^2`.
4. **CHECK 1's claimed sensitivity is overstated.** “any error in the eigenstrain factor ... or ... sign/factor conventions ... would make the integral differ” (`submission/ANSWER.md:104–108`). The code integrates the square of an elastic stress trace; it contains no eigenstrain coefficient and cannot test its sign or prefactor. Its successful integral cannot validate the short-time application of reciprocity.
5. **CHECK 2's physical validation is overstated.** “independent 2-D simulation of the short-time diffusion” and “the wrong diffusivity ... would shift C0” (`submission/ANSWER.md:110–132`). The calculation solves scalar heat flow with diffusivity fixed to one, without equilibrium, constitutive coupling, or recovery of pressure from fluid content. It neither tests the supplied physical diffusivity nor validates the Biot crack-volume coefficient. Its fitted leading coefficient is valid evidence for its own surrogate only.
6. **Slit and reference crack are misaligned in CHECK 2.** `submission/solve.py:62–71` imposes the slit on row `n//2`; the actual coordinate is `y = −0.0005`, while the reference crack is at zero. The comment “grid row at ... y = 0” is inaccurate. This is a small numerical geometry error in the check, not an additional endpoint error.

**Presentation: 1.** “an incompressible-elastic-type response” (`submission/ANSWER.md:14`) is imprecise: undrained elasticity is incompressible only when `nu_u = 1/2`, whereas the example has `nu_u = 0.4`. The equations written alongside it use the correct undrained ratio.

## 5. Checks

There are **four reported checks**: two numerical checks and two expressly claimed analytical verifications. **Two have independently computed outcomes that could disagree with their predictions.** No check determines the drained-tip intensity or the actual coupled short-time coefficient. Neither numerical check asserts a tolerance; disagreement would have to be noticed in its printed output.

| Reported check | Tests the deciding quantity? | Could fail? |
|---|---|---|
| CHECK 1: numerical `int sigma_kk^2 dA` versus `2*pi*L^2` | Tests an elastic-reference norm and supports the correct endpoint-volume difference; does not test the disputed tip or short-time coefficient, or reciprocity prefactors | Yes: an incorrect reference field or numerical integral can disagree |
| CHECK 2: scalar slit diffusion, weighted integral, fitted coefficient versus 16 | Tests the face-layer coefficient of a surrogate; omits the deciding poroelastic coupling and leading tip term | Yes for the surrogate: a wrong diffusion computation or predicted surrogate coefficient can disagree; it can pass while the physical answer is wrong |
| Analytical `int_0^inf erfc(s) ds = 1/sqrt(pi)` (`ANSWER.md:134`) | A correct auxiliary identity, not a check of the solution | No independent verification is shown; the quoted identity is unchanged by solution errors |
| Analytical consistency of `c` and the `zeta` equation (`ANSWER.md:135–136`) | Correct elimination, but does not verify that pressure obeys the surrogate heat equation | No separate comparison is shown; it restates the algebra used to obtain the supplied coefficient |

The automatic volume-ratio and Griffith-limit checks also pass. They were generated by the grading system and are not added to the submission's check count.

## 6. Claims and reproduction

**Supported as numerical reproduction.** `AUTO_CHECKS.json` reports successful execution and all nine claims printed. A fresh run of `python3 solve.py` in `scratch/submission/` exited zero and reproduced the complete original `results.json`, including formula strings. Its two numerical check outputs match the pasted output, including 6.277848, relative error `8.494e-04`, 39 fit samples, `C0 = 16.1118`, and RMS `4.84e-05`. There are **no numerical or pasted-output mismatches**. Evidence is in `scratch/submission/rerun_output.txt` and `scratch/audit_results.json`.

The scientific accuracy and limitations of the checks are separately assessed above. The two analytical verification claims are prose-only; the reproduction command implements the two numerical checks it promises.

## 7. Questions and missing information

Four clarification/assumption bullets are listed under “Questions / assumptions” (`submission/ANSWER.md:140–155`); there are no outstanding information requests. Reading `sigma_n` as total traction is sensible and agrees with the statement, but does not establish the claimed zero loading-2 long-time response. The alternative cavity loading is correctly identified as the sum of the two loadings; its quoted endpoint predictions inherit the errors/challenge already recorded. Distinguishing the initial interior pressure from the imposed face value is sensible. Unit loading and plane strain merely repeat supplied information. “Nothing missing” is stated at lines 159–161.

## 8. Source

**No source paper or author is identified for this poroelastic problem.** The phrases “the Westergaard solution” (`submission/ANSWER.md:25`) and “Betti reciprocal theorem” (line 21) name an elastic field and a theorem, but do not identify the problem's source paper or its authors. No Detournay/Cheng citation or other poroelastic source attribution is delivered. Thus `source_named` is `none`.

## 9. Challenges and assessment of the key

**One substantive, checkable challenge:** under literal zero total traction, loading 2 has a uniform-pressure, stress-free long-time limit, so `K2_inf = V2_inf = 0`. The argument is: “in 2-D there is no bounded steady state” with the given pressure conditions; locally `p -> 1`; then “sigma = 0 ... [satisfies] all boundary conditions (unique)” and uniform free strain cannot open a slit (`submission/ANSWER.md:65–72,140–147`). This is an explicit limiting argument, so its affected answers are marked contested and its two disputed endpoint values are not counted as errors.

**Checkability and evidence:** one can retain the far field in finite-radius or increasing-time coupled solutions and measure the local limiting stress and crack response. The absence of a bounded global steady pressure field is correct. The asserted uniqueness of the stress-free local limit does not follow: taking the local pressure limit before retaining far-field mechanical matching loses a constraint.

A concrete counterexample to that inference is computed in `scratch/audit.py`. Let pressure be one inside radius `R` and zero outside in an uncracked plane-strain solid. With `lambda = 2G*nu/(1−2nu)`, the matched displacement is `u = a*r` inside and `u = a*R^2/r` outside, where `a = alpha/[2(lambda+2G)]`. Displacement and radial traction match at `R`, and the fields decay at infinity for each finite `R`. As `R -> infinity` at fixed position, pressure tends to one but both in-plane total stresses remain `−eta`, where `eta = alpha(1−2nu)/[2(1−nu)] = 0.33375`. Thus locally uniform pressure does not imply zero stress. Cancelling the background normal stress on a traction-free slit requires the closing elastic correction, giving the key's `K2_inf = −eta` and `V2_inf = −eta*pi*(1−nu) = −0.8388052385` in the stated units.

This pressure-patch calculation tests the submission's limiting inference; it is not an independent transient solution of the original slit problem. It supplies evidence favoring the key's retained confining stress. Against the proposed stress-free end-state assumption, (b) is conditionally right and (c) remains partly right; against the stated problem and key, both are partly right. I find no demonstrated error in the corrected key. The separate loading-1 order-of-limits ambiguity is already acknowledged by the key and earns the specified partial credit.

```json
{"answers": {"a": "partly", "b": "contested", "c": "contested"},
 "supporting_state": "right",
 "method": "elastic limits plus uncoupled pressure diffusion and Betti reciprocity", "method_class": "general_or_wrong",
 "errors": {"decisive": 3, "side": 6, "presentation": 1},
 "checks_reported": 4, "checks_could_fail": 2,
 "claims_supported": "yes",
 "questions_listed": 4,
 "source_named": "none",
 "challenges": [{"claim": "Literal traction-free loading 2 has a stress-free uniform-pressure long-time limit, giving K2_inf = V2_inf = 0", "checkable": true}]}
```
