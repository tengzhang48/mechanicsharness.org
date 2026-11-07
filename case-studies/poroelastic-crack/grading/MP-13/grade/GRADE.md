# Grading profile

## 1. Answers

Grades use the **P374 corrected key**, with its 1% nonzero-value bands and absolute tolerance 0.01 for `K2_0`.

- **(a) Partly right.** Both endpoint volumes and the dimensional long-time K formula are correct; the claimed constant K misses the drained tip limit, and the short-time growth coefficient is wrong.
- **(b) Wrong.** All three formulas give positive opening/intensity; the key requires `K2_0 = 0`, negative long-time K, and negative long-time volume (closure).
- **(c) Partly right.** Only two of the nine delivered numerical fields satisfy the stated bands; an additional normalization error affects every K field.

| Field | Submitted JSON | Corrected key | Within band? |
|---|---:|---:|---|
| K1_0 | 1.772454 | 0.807593 | No |
| K1_inf | 1.772454 | 1.000000 | No |
| V1_0 | 1.884956 | 1.884956 | Yes |
| V1_inf | 2.513274 | 2.513274 | Yes |
| V1_short_coeff | 1.128379 | 2.440220 | No |
| p_faces_0 | −0.118667 | +0.749064 | No |
| K2_0 | 1.183113 | 0.000000 | No |
| K2_inf | 1.577484 | −0.333750 | No |
| V2_inf | 2.236814 | −0.838805 | No |

Dividing the submitted K values by √π at L = 1 gives `(1, 1, 0.6675, 0.89)`. This repairs the long-time loading-1 number only. The initial loading-1 formula, K = √(πL), deserves the key's partial credit as the undrained **outer** intensity under the reverse order of limits. The submission never distinguishes that limit from the drained tip intensity, whose corrected normalized value is 0.807593. Its assertion that K is constant for all time remains wrong. The paper's superseded 0.75 and 1.354055 values are not used as the grading targets.

## 2. Supporting state

**Wrong for the requested induced face-adjacent pressure.** The undrained constraint gives

`p = −(ν_u − ν) σ_kk / [α(1 − 2ν)]`.

The loading-1 undrained outer face field has σ_kk = −2, giving `p_faces_0 = 2(ν_u − ν)/[α(1 − 2ν)] = +0.7490636704`. The submitted −0.1186666667 has the wrong sign, magnitude, and parameter dependence. Its formula even retains nonzero pressure as ν_u approaches ν.

Correct supporting pieces include unchanged shear modulus G, the bulk undrained compliance, the two loading-1 volume formulas, their ratio 4/3, and the diffusivity expression in `solve.py:44–47`. These do not supply the missing tip analysis.

## 3. Method

**Elastic endpoints with assumed relaxation and effective traction — general_or_wrong.** The solver uses elastic crack formulas, assumes a common scalar time factor for both loadings, substitutes effective tension αp for loading 2, and assigns an unsupported tip-layer volume coefficient. It supplies no transient solution.

Elastic bulk endpoints are useful, but linearity does not make the two distinct boundary loadings share a relaxation function. The method omits the drained tip limit, the pressure-gradient contribution to equilibrium, and the diffusion/reciprocity calculation needed for short-time growth, including the corrected tip contribution. The discarded finite-difference idea in the process log is not the delivered method.

## 4. Errors

Counts below refer to distinct incorrect formulas or assertions; repetitions across files are counted once. Related causes and consequences are stated together where appropriate. Totals: **9 decisive, 5 side, 2 presentation**.

| ID | Class | Quote or location | Error and consequence |
|---|---|---|---|
| D1 | decisive | `solve.py:51,99–106`; `results.json:2–3,8–9` | K is computed with `S = sqrt(pi*L)` and written without division by the requested unit √(πL). All four K fields have the wrong numerical normalization. This changes answer (c), so it is decisive rather than merely a unit label error. |
| D2 | decisive | `ANSWER.md:47–49`: “time-independent” | Stiffness independence of an elastic traction crack is extended to the evolving poroelastic tip field. K1_0 misses the drained tip intensity and its limit order; the correct long-time formula does not justify constant K. |
| D3 | decisive | `solve.py:61–63`; `ANSWER.md:53–56` | `A = (8/(3√π))(1−ν_u)/(1−ν)L²/G` has no demonstrated diffusion/opening integral and gives 1.128379 rather than 2.440220. A remains 1.504506 in normalized units as ν_u approaches ν, when the drainage-induced volume growth should disappear. The asserted tip-only relaxation also omits the face contribution. |
| D4 | decisive | `solve.py:65–69`; `ANSWER.md:105–109` | The final induced-pressure formula is negative, lacks ν_u, and scales with α instead of the required inverse α. It contradicts the undrained constraint and the opening-load face field. |
| D5 | decisive | `ANSWER.md:81`; `solve.py:67` | The purported checked relation `σ'_kk = +α(1−2ν_u)p/(ν_u−ν)` has the wrong sign. From ζ = 0 and σ'_kk = σ_kk + 2αp, the result is **−**α(1−2ν_u)p/(ν_u−ν). It is cited as the basis of the incorrect pressure calculation, and its claimed successful check is false. |
| D6 | decisive | `ANSWER.md:36–43,64`; `solve.py:75` | The common-Φ inference from linearity produces nonzero `K2_0 = α√(πL)(1−ν_u)/(1−ν)`. The hydraulic boundary loading has zero initial crack-face displacement and zero initial K under the key's limit. Uniform boundary data alone do not make the two responses proportional. |
| D7 | decisive | `ANSWER.md:60–63`; `solve.py:74` | Treating αp as a complete opening traction gives `K2_inf = +α√(πL)`. Pressure acts throughout the poroelastic field, producing the back stress; the correct result is `−α(1−2ν)√(πL)/[2(1−ν)]`. |
| D8 | decisive | `ANSWER.md:65`; `solve.py:76` | `V2_inf = +πα(1−ν)L²/G` predicts opening instead of closure. The correct formula is `−απ(1−2ν)L²/(2G)`. |
| D9 | decisive | `ANSWER.md:26–27,61,119–120` | The drained endpoint is described as `p ≡ 0 away from the crack`, and loading-2 pressure as draining away. Loading 2 maintains p = 1 on the faces for every t > 0; its long-time pressure-induced stress cannot be discarded. This incorrect field assumption underlies its elastic opening substitution. |
| S1 | side | `ANSWER.md:23–24`: `σ'_ij,j = −α p,i` | With σ'_ij = σ_ij + αpδ_ij and total-stress equilibrium, the derivative is **+**αp,i. The delivered code never solves this equilibrium equation, so this additional sign error does not alter its implemented numbers. |
| S2 | side | `ANSWER.md:93–100` | The proposed decomposition contradicts its own full hydraulic-loading anchor. In normalized units its final formulas sum to K(∞) = 1 + 0.89 = **1.89**, while the cited anchor says 1; initially they sum to **1.6675**, while the anchor says 0.75. No consistent superposition establishes the claimed split. |
| S3 | side | `ANSWER.md:118–122`: “K decreases with time” and “opposite time trends” | Its own loading-2 formulas increase from 0.6675 to 0.89 after normalization, and its loading-1 K is constant. The qualitative conclusion contradicts its formulas and is not a check of them. |
| S4 | side | `submission/PROCESS.md:3`: `drained sigma'_kk=-alpha/(1-nu)p` | There is no such universal drained algebraic constraint: steady drainage does not impose ζ = 0, and a mechanically loaded drained crack can have nonzero stress trace with p = 0. The relation is not used by the final code. |
| S5 | side | `submission/PROCESS.md:3`: “crack modulus ratio ... =9/4” | That expression is the ratio of longitudinal constrained moduli, not crack stiffnesses. The endpoint crack-volume formulas imply a crack-stiffness ratio `(1−ν)/(1−ν_u) = 4/3`. The 9/4 claim is unused in the final calculation. |
| P1 | presentation | `ANSWER.md:109`: “closing traction”; `solve.py:54`: `q = -1` | For the elastic opening formula quoted by the solver, the positive opening load is q = −σ_n = +1. Calling loading 1 a closing traction, and assigning q = −1 in that convention, conflicts with its correctly positive loading-1 volumes. The numerical pressure sign error is separately counted in D4. |
| P2 | presentation | `ANSWER.md:3`: “prints everything below” | The executable prints nine numbers and a file-write message; it prints no formulas, reasoning, questions, or symbolic checks. |

The first-pass volume factor 2, spurious shear-modulus/R(0) bookkeeping, and provisional short-time rate in `submission/PROCESS.md:2` are explicitly marked wrong and withdrawn; `ANSWER.md:76–80` acknowledges the corrections. They are recorded here as corrected history and receive no additional final-error counts.

## 5. Checks

There are **five numbered reported checks** (`ANSWER.md:72–89`). The process log repeats the symbolic claims and adds intermediate assertions covered above, without separate test procedures. The automatic grading checks and this grader's additional limit checks are not counted as solver-reported checks.

| Reported check | Deciding quantity tested? | Could fail for an incorrect result? | Assessment |
|---|---|---|---|
| 1. Undrained compliance | Supports bulk opening/endpoint volume; does not test tip K | Yes | Direct elimination independently confirms Hooke compliance with ν_u and unchanged G. No executable SymPy check is delivered. |
| 2. Ellipse volume integral | Yes, the volume factor for the assumed opening | Yes | Independent symbolic integration confirms πq(1−ν)L²/G. It does not validate the hydraulic loading or the short-time correction. |
| 3. Undrained effective-stress trace | Yes, pressure/sign | Yes; it should have failed | The reported positive-sign identity is false. Independent elimination gives the negative sign shown in D5. |
| 4. Diffusivity identity | No endpoint, tip, or short-time coefficient is tested | No independent failure mechanism shown | The code copies the statement's c verbatim; neither the claimed deduction nor an independent comparison is delivered. |
| 5. Run/output check | Reproducibility only | No for a wrong mechanical solution | The run succeeds and reproduces the files, including their incorrect values. A missing NumPy dependency would test the environment, not the solution. |

**checks_could_fail = 3**, counting the independently falsifiable mathematical targets in checks 1–3. This is not a count of three valid passed checks: check 3's claimed pass is incorrect. The executable contains no assertions or symbolic validations.

## 6. Claims

**Partly supported overall.** `AUTO_CHECKS.json` reports exit 0 and all nine numerical claims printed. A fresh run of `python3 solve.py` in `scratch/rerun/` confirms exit 0, exact equality of regenerated and delivered JSON, and agreement of all printed values to their displayed precision. The original submission was left unchanged. `scratch/verify.py` and `scratch/rerun_stdout.txt` preserve the local verification and output.

Numerical reproducibility is complete. The delivered program does not perform the claimed symbolic checks, and the claimed effective-stress trace identity is false. The stated K units conflict with the computed fields, and the prose's “1.00 at both times” (`ANSWER.md:49`) conflicts with the table/JSON value 1.7724538509. The AUTO_CHECKS failures and volume-ratio success are all reproduced.

## 7. Questions

**Four question/uncertainty bullets** occur in `ANSWER.md:93–114`.

- **Loading split/coupling:** a relevant issue to investigate, but the given total-traction and pressure boundary conditions determine the coupling. The α versus α/[2(1−ν)] alternatives are unsupported assumptions, not missing input; neither gives the key's signed back-stress coupling.
- **Induced-pressure bookkeeping:** the location before diffusion is a sensible clarification. The claimed freedom to obtain different pressures by bookkeeping is unsupported; the constitutive constraint fixes the relation and sign.
- **Short-time prefactor:** sensible uncertainty and a proposed numerical verification, but no BEM/FEM calculation or coefficient extraction is delivered.
- **ν_u > 1/2:** already answered by the stated range ν < ν_u ≤ 1/2, and irrelevant to the supplied value 0.4.

## 8. Source

The solver names **“Rice & Cleary, JGR 81, 409, 1976”** (`solve.py:9`), also “Rice & Cleary (1976)” (`ANSWER.md:36–37`). It does not name the key's Detournay–Cheng paper. No supplied source passage establishes the claimed common relaxation factor or 8/(3√π) coefficient; no external lookup was performed.

## 9. Challenges

Two statements amount to limited ambiguity/incompleteness challenges. Neither establishes a well-argued challenge qualifying for a contested grade.

1. **Claim:** the loading decomposition is “not unique from R&C alone,” and the intended effective coupling might be α/[2(1−ν)] rather than α (`ANSWER.md:93–104`). **Argument:** a shared time factor is assumed from uniformity/linearity, followed by two explicit candidate couplings. **Checkable: yes**, through the supplied formulas, superposition, and boundary-equilibrium checks. In particular, adding its own two endpoint K formulas contradicts its stated hydraulic-loading anchor, as S2 demonstrates. Linearity supplies superposition, not the assumed common response function. Against the key, (b) is wrong; even using its proposed alternative coupling leaves nonzero initial K and positive long-time K/volume, so it does not rescue (b).

2. **Claim:** a different, “equally plausible” pressure/body-force bookkeeping changes `p_faces_0` by an O(1) factor (`ANSWER.md:105–109`). **Argument:** only the assertion of another bookkeeping choice; no alternative pressure field or formula is supplied. **Checkable: no as an ambiguity claim**, because there is no second construction to reproduce. Its actual reported formula is checkable and fails the undrained constraint. Against the key the supporting pressure is wrong; there is no defined alternative solution to grade against this challenge.

The initial K limit ambiguity is acknowledged by the key and earns the partial credit stated in section 1; the solver does not articulate that limit-order challenge. No additional defect in the key was established by these local checks.

```json
{"answers": {"a": "partly", "b": "wrong", "c": "partly"},
 "supporting_state": "wrong",
 "method": "Elastic endpoints with assumed relaxation and effective traction", "method_class": "general_or_wrong",
 "errors": {"decisive": 9, "side": 5, "presentation": 2},
 "checks_reported": 5, "checks_could_fail": 3,
 "claims_supported": "partly",
 "questions_listed": 4,
 "source_named": "Rice & Cleary, JGR 81, 409, 1976",
 "challenges": [{"claim": "Loading decomposition/coupling is not determined by the cited Rice–Cleary anchor", "checkable": true}, {"claim": "Equally plausible pressure bookkeeping changes p_faces_0 by an O(1) factor", "checkable": false}]}
```
