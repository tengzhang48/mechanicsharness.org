# Grading profile

## 1. Answers

Grading uses the **P374 corrected key**, not its retained historical formulas: `K1_0 = 0.807593` and `V1_short_coeff = 2.440220`. The bands are 1% for nonzero values and 0.01 for `K2_0`.

| Question | Judgment | Reason |
|---|---|---|
| (a) | **Partly right** | Both endpoint volumes and the long-time K are right; initial K uses the undrained outer limit rather than the drained-tip limit, and the leading volume-growth coefficient is wrong. |
| (b) | **Right** | All three requested limits, including the negative signs of the long-time K and volume, agree with the key; there is nevertheless an error in the supporting steady pressure field. |
| (c) | **Partly right** | Seven of nine numbers satisfy the bands; `K1_0` receives the key's explicit partial credit for the alternative order of limits, while `V1_short_coeff` is wrong. |

| Quantity | Submitted | Corrected key | Assessment |
|---|---:|---:|---|
| `K1_0` | 1 | 0.807593 | Partly: valid outer limit, 23.82% above the stipulated tip limit |
| `K1_inf` | 1 | 1 | Right |
| `V1_0` | 1.884955592 | 1.884955592 | Right |
| `V1_inf` | 2.513274123 | 2.513274123 | Right |
| `V1_short_coeff` | 1.128379167 | 2.440220 | Wrong: 53.76% low |
| `p_faces_0` | 0.749063670 | 0.749063670 | Right |
| `K2_0` | 0 | 0 | Right |
| `K2_inf` | -0.33375 | -0.33375 | Right |
| `V2_inf` | -0.838805239 | -0.838805239 | Right |

## 2. Supporting state

**Right for the requested supporting quantities.** The pre-diffusion adjacent-solid pressure `p_faces_0 = 2(nu_u-nu)/[alpha(1-2nu)]` is correct. The undrained relation `p = -(nu_u-nu)s/[alpha(1-2nu)]`, face trace `s = -2`, elliptical endpoint openings, and coupling constants `gamma = 0.6675`, `eta = gamma/2 = 0.33375` are also correct. The distinction between this adjacent-solid pressure and the imposed zero face pressure is stated correctly (`ANSWER.md:130`). The additional simulated transient and proposed steady pressure fields have errors listed below; they are not requested intermediate outputs.

## 3. Method

The solver combines elastic Sneddon/Westergaard limits, an eigenstrain decomposition, a local one-dimensional erf drainage layer, and Chebyshev crack calculations with finite-difference diffusion and FFT eigenstress checks. Elastic endpoint calculations are appropriate. The local layer is insufficient to determine the whole-crack volume coefficient: reciprocal/global matching and the leading tip contribution are missing. The undrained elastic K is substituted for the drained-tip limit. The numerical transient is both unnecessarily general for the requested asymptotics and incorrectly implemented. Overall: **general_or_wrong**, with correct components rather than a valid complete alternative.

## 4. Errors

Counts refer to distinct errors in the delivered solution. Repeated appearances are grouped. Historical mistakes explicitly repaired before delivery are not counted again. Errors confined to extra verification code are side errors because `solve.py` evaluates the reported answers separately.

### Decisive — 3

1. **Wrong initial-tip K and unwarranted time independence.** “SIF ... is fixed by tractions alone” (`ANSWER.md:29`; `solve.py:21`, `solve.py:55`). This is valid for the homogeneous elastic outer solution, but pore-pressure-induced stresses change the poroelastic tip intensity. The key specifies 0.807593 for the tip-first limit. The value 1 merits partial credit under the reverse order; the solver neither distinguishes these limits nor resolves the drained-tip zone. The stronger claim `K1(t) = ... for all t` is also unjustified.
2. **Incorrect conversion of the local drainage layer into whole-crack volume.** “gives uniform opening change ... hence ΔV = 2LΔw” (`ANSWER.md:35`; `solve.py:32`). The locally constrained strain integral is treated as the complete global opening response without reciprocal/global matching. Even the key's face contribution alone is `16(1-nu_u)(nu_u-nu)/[sqrt(pi)(1-nu)] = 1.354055`, rather than 1.128379.
3. **Leading tip contribution incorrectly discarded.** “tip-region corrections are O(ct/L²)” (`ANSWER.md:37`, repeated at `ANSWER.md:132` and `results.json:15`). The corrected key has an additional **1.086165 at order sqrt(ct/L²)**. This omission is distinct from the incorrect face coefficient: fixing either one alone leaves the total wrong.

### Side — 10

1. **False negative-diffusivity assertion.** “The alternative reading (σ_kk = 3D trace) yields ... (negative) diffusivity” (`ANSWER.md:11`). No calculation is supplied. Even literally extending these coefficients to three dimensions and imposing `e_zz = 0` gives `c_3 = 2*kappa*G*(nu_u-nu)/[alpha²(1-2nu)(1-2nu_u)]`, which is positive for `nu_u < 1/2` and equals 4.208223 at the given parameters with `kappa = G = 1`. It differs from the stipulated c; it is not negative. The intended indices 1,2 are already explicit.
2. **Wrong loading-2 steady pressure boundary problem.** “p_ss is harmonic (strip potential, 1 on faces, →0 at ∞)” (`ANSWER.md:42`; `verify.py:96`). This strip potential imposes `p = 0` on the intact plane outside the crack. Symmetry and pressure continuity require `p_y = 0` there. Its derivative at `x = 2L, y = 0+` is `2/(3*pi*L)`, not zero. Thus it is not the crack's pressure solution, although inserting its uniform face trace produces the right reported loading-2 limits.
3. **Incomplete transient governing equation.** The transient evolves `p_t = c*laplacian(p)` (`partC.py:2`, `partC.py:95`, `partC.py:151`). Define `H = s + gamma*p`, `beta = alpha(1-2nu)`, and `D = alpha²(1-2nu)²(1-nu_u)/[(nu_u-nu)(1-nu)]`. Mass balance gives `p_t = c*laplacian(p) - (beta/D)*H_t`. Spatial harmonicity of H does not imply `H_t = 0`. The answer acknowledges this approximation (`ANSWER.md:116`) but never establishes it for the crack transient.
4. **Neumann rows are not enforced.** Boundary rows are inserted into `Lap`, then transformed by `Aimp = I-dt*Lap` and given zero right-hand sides (`partC.py:72`, `partC.py:77`, `partC.py:95`, `partC.py:150`). A intended row `p_0-p_1 = 0` instead becomes `(1-dt)*p_0 + dt*p_1 = 0`. At `dt = 5e-6`, neighbor pressure 1 produces boundary pressure -0.000005000025, rather than 1. Both the intact-plane and x-symmetry conditions are wrong.
5. **Reflection is shifted by one cell.** The arrays are reversed onto a grid containing zero (`partC.py:128`, `partC.py:129`). The point `x = -h` receives the value at `x = 0`, rather than `x = +h`; likewise in y. For the specified FFT grid, `h = 0.00390625`. The resulting field is reflected about `-h/2`, rather than the stated axes.
6. **Wrong crack-volume integration for nonuniform loads.** `V = -(pi*L/2)*sum(a[odd]/n[odd])` (`verify.py:38`; `partC.py:45`). From the code's own opening `w = -sum(a_n*sin(n*theta)/n)`, orthogonality gives **`V = -(pi*L/2)*a_1`**. A manufactured load `q = 1 + 0.5*U_2(x/L)` gives coded V = 2.932153, but integration of that same opening gives 2.513274. Uniform-load checks cannot expose this error. It contaminates the extra transient V and fitted A.
7. **Correction-field K reported as total K.** The volume includes its baseline elastic contribution, while the transient K does not (`partC.py:143`, `partC.py:147`, `partC.py:157`). The saved and rerun values are `Kc = -0.18009` initially and -0.18629 to -0.15871 later. Adding the elastic baseline 1 gives the positive numbers paraphrased in `ANSWER.md:115` and `ANSWER.md:137`; the delivered output and NPZ do not contain those total values. Their variation also does not validate a constant K.
8. **Missing L in the stress-intensity quadrature.** The substitution `x = L*cos(theta)` requires `dx = -L*sin(theta)*dtheta`, but that factor L is absent (`verify.py:42`; `partC.py:47`). With `L = 2`, uniform unit loading gives normalized K = 0.5 instead of 1. This is hidden by every delivered run using `L = 1`.
9. **Claimed face-stress verification is absent.** “numerically confirmed” and “it also confirms the face values σ_xx = σ_yy = -σ₀” (`ANSWER.md:38`, `ANSWER.md:78`). The code computes only opening, V, K, and **σ_yy ahead of the tip**, at `x > L` (`verify.py:80`). It computes neither face σ_xx nor face σ_yy. The face values themselves are correct, but this stated evidence is not delivered.
10. **Incorrect attribution of the steady FFT discrepancy to tip aliasing.** “underestimates ... due to tip-jump aliasing at N=2048” (`ANSWER.md:103`). At fixed box size 4, the sample trace converges from -0.256427 at N=256 to -0.256423 at N=2048, retaining a 23.17% error. Enlarging the box reduces the error to 11.55% and 5.77%. This demonstrates a substantial finite-box/periodic-field error; resolution of the tip jump does not explain the quoted discrepancy. Likewise, the transient discrepancy cannot be assigned just to aliasing when its boundary conditions and volume operator are wrong.

### Presentation — 2

1. **Stale opposite-sign comment:** `verify.py:116` says `q = -sigma_yy^sm = gamma/2`, while the correct executed input is `q = -gamma/2` at line 117.
2. **Stale script description:** `verify.py:9` describes a Part C transient implemented by that script, but the delivered file ends after Part B. The transient is in `partC.py`, correctly named in the answer.

## 5. Checks

Counting convention: **seven reported checks**, treating the three explicitly separate checks within numbered group 2 as separate checks. Multiple observables from one run are not counted again. “Could fail” describes the stated comparison's ability to detect an error in its tested component, not successful validation of the whole solution.

| Reported check | Deciding quantity and coverage | Could fail? |
|---|---|---|
| 1. Uniform crack versus closed forms | Tests elastic w, V, K, and ahead-tip stress; rerun matches. Does not test poroelastic tip drainage, nonuniform-volume integration, or the claimed face stresses. | Yes |
| 2(i). Harmonic-mode FFT versus analytic stress | Tests the smooth-field kernel relevant to loading-2 stress/sign. Arrays are shown, but no code or exact grid/mode inputs are delivered for reproducing this run. | Yes, as a reported numerical comparison |
| 2(ii). Disk inclusion | Tests the kernel's sign and magnitude; rerun gives -0.32965 versus -0.33375, a 1.23% discrepancy. Does not solve the crack hydraulic boundary problem. | Yes |
| 2(iii). x-uniform slab limit | The independent 1D equilibrium limit requires σ_yy = 0; the stated kernel satisfies it at `kx = 0`. This is an analytical check, with no separate run delivered. | Yes, for an incorrectly oriented kernel |
| 3. Loading-2 steady crack calculation | Uniform imposed `q = -gamma/2` gives the reported K and V. It tests the elastic conversion, while assuming the decisive back-stress. The separate FFT-derived route fails: K = -0.250706 and V = -0.643652. | Yes, with limited coverage |
| 4. Loading-1 transient | Attempts the deciding initial K and growth coefficient. It actually fails: A_fit = 2.69359 versus 1.12838, intercept 2.00182 versus 1.88496, and incorrect K reporting. The governing and measurement bugs prevent quantitative validation. A 0.7% endpoint V agreement does not validate the pipeline. | Yes; failure is disclosed |
| 5. Diffusivity consistency | Conditional constitutive/compatibility algebra is correct: `c = 2*kappa*G/D`. It can discriminate an inconsistent trace/coefficient convention. The submitted script stores this as a textual claim, rather than executing a comparison. It does not establish `H_t = 0`. | Yes, as an analytical consistency check |

**Checks that could fail: 7.** None provides a successful independent validation of the required drained-tip K or complete leading volume coefficient.

## 6. Claims and reproduction

**Partly supported overall.** `AUTO_CHECKS.json` reports successful reproduction of all nine numerical claims. Rerunning `python3 solve.py` in `scratch/submission/` confirmed all nine printed values and an exactly identical parsed `results.json`. Reproducibility does not make the two incorrect results correct.

Rerunning `verify.py` reproduces the printed elastic and imposed-back-stress comparisons, including the failed raw-FFT route. Rerunning `partC.py` reproduces the submitted transient output and NPZ arrays, apart from timing text. The quoted V values and A_fit = 2.69359 are supported. The written positive K range is not what the delivered script prints or saves; the baseline is omitted as described above. Harmonic-mode arrays and a separate executed slab check cannot be reproduced from the delivered scripts. The diffusivity check is supplied as reasoning/text, not runnable verification. The attribution of failed comparisons to aliasing is not supported by the box/resolution diagnostic.

Reproduction logs and independent diagnostics are retained under `scratch/`: `submission/solve_rerun.txt`, `submission/verify_rerun.txt`, `submission/partC_rerun.txt`, `audit.py`, and `audit_results.json`.

## 7. Questions and missing information

**Six numbered items** appear under “Questions / assumptions” (`ANSWER.md:127`), although they are statements rather than questions. The loading normalization and adjacent-solid versus face-pressure distinction are sensible and already specified. The trace convention is explicit in the problem, so no clarification was needed. The claimed order of tip corrections and constant-K inference are errors, not missing information. The contact caveat is sensible for interpreting negative linear openings; contact is outside the supplied model. The later uncertainty notes correctly disclose the unresolved A discrepancy and incomplete transient sampling. No necessary material parameter or boundary condition is missing, and no actual order-of-limits question is raised.

## 8. Source

The submission names authors for its model and classical checks: **“Rice–Cleary-type plane-strain poroelasticity in 2D form”** (`ANSWER.md:6`), “Sneddon/Westergaard closed forms” (`ANSWER.md:23`), and “Eshelby disk inclusion” (`ANSWER.md:86`). Credit these attributions. It does **not** identify a source paper for this particular problem or name Detournay and Cheng, the key's source authors.

## 9. Challenges and key issues

The solver's interpretation challenge is **“The alternative reading (σ_kk = 3D trace) yields an inconsistent (negative) diffusivity”** (`ANSWER.md:11`). Its argument is that matching c selects the in-plane trace. This is checkable by eliminating σ_zz with `e_zz = 0` and recomputing the storage coefficient, but no supporting 3D calculation is supplied and the negative-sign claim fails that calculation. It is therefore an error, not a well-argued challenge qualifying for a contested grade. The contact statement is a conditional caveat about adding different physics, rather than an argument that the supplied linear problem is wrong. No published result or author is challenged.

For the acknowledged ambiguity of initial K, the two grades are explicit: **against the key's tip-first limit, K1_0 is not right; against the reverse outer-limit reading, K1_0 = 1 is right.** The key expressly assigns partial credit, which is used here. The solver supplies no valid challenge to the growth coefficient; (a) and (c) remain partly right under either initial-K reading.

One minor transcription issue in the key is independent of the submission: its physical opening–K relation at `KEY.md:61` prints `sqrt(pi*G)` where the prefactor should be `G*sqrt(pi)`. From `w(r) = [sqrt(8)*(1-nu)/G]*K*sqrt(r/pi)`, inversion gives `K = G*sqrt(pi)/[sqrt(8)*(1-nu)]*lim(w/sqrt(r))`. The printed form is dimensionally inconsistent for dimensional opening and tip distance. Correcting this factor does not change the supplied normalized benchmark values or any of the judgments above. Grading with either transcription therefore gives the same (a)/(b)/(c) profile; the historical 0.75 and face-only 1.354055 remain superseded by P374.

```json
{"answers": {"a": "partly", "b": "right", "c": "partly"},
 "supporting_state": "right",
 "method": "Elastic limits, local 1D drainage, FD/FFT and Chebyshev verification", "method_class": "general_or_wrong",
 "errors": {"decisive": 3, "side": 10, "presentation": 2},
 "checks_reported": 7, "checks_could_fail": 7,
 "claims_supported": "partly",
 "questions_listed": 6,
 "source_named": "Rice–Cleary-type plane-strain poroelasticity in 2D form",
 "challenges": [{"claim": "A 3D-trace interpretation gives negative diffusivity", "checkable": true}]}
```
