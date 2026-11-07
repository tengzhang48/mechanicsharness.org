# Grading profile

## 1. Answers

The corrected P374 values in KEY.md take precedence over its retained paper formulas: `K1_0 = 0.807593` and `V1_short_coeff = 2.440220`. Nonzero reference values have a 1% tolerance; `K2_0` has an absolute tolerance of 0.01.

- **(a) Partly right:** the two volumes and long-time intensity are right, but the initial intensity is wrong and the leading small-time growth is not obtained.
- **(b) Wrong:** both requested intensities and the long-time volume are wrong; the volume has the wrong sign as well as the wrong magnitude.
- **(c) Partly right:** three of the nine requested numerical quantities pass; six fail, as AUTO_CHECKS.json reports.

All entries below use the problem's normalized units.

| Quantity | Submission | Corrected key | Verdict |
|---|---:|---:|---|
| `K1_0` | 1.25 | 0.807593 | Wrong |
| `K1_inf` | 1 | 1 | Right |
| `V1_0` | 1.8849555922 | 1.8849555922 | Right |
| `V1_inf` | 2.5132741229 | 2.5132741229 | Right |
| `V1_short_coeff` | 0, explicitly a placeholder | 2.440220 | Wrong / missing |
| `p_faces_0` | 2 | 0.7490636704 | Wrong |
| `K2_0` | -1.25 | 0 | Wrong |
| `K2_inf` | -1 | -0.33375 | Wrong |
| `V2_inf` | 2.5132741229 | -0.8388052385 | Wrong |

The ambiguity acknowledged by the key does not rescue `K1_0`: the reverse order of time and tip limits gives 1, also different from 1.25. The submission neither identifies that order-of-limits issue nor offers a supported alternative to the corrected key.

## 2. Supporting state

**Wrong for the requested supporting state.** The separately requested intermediate quantity, `p_faces_0`, is wrong. The Sneddon opening `w = 2(1-nu)*sqrt(L^2-x^2)/G` and its undrained counterpart are correct (ANSWER.md:42–45, solve.py:103–106), and explain the volume credit in (a) and (c). The supplied diffusivity expression is also correct. The additional statement `V2_0 = 0` is consistent with the key, although its stated delta-response explanation is wrong.

The requested induced pressure is wrong. For the undrained outer state, setting fluid-content variation to zero in the given constitutive relation gives `p = -(nu_u-nu)*sigma_kk/[alpha*(1-2*nu)]`; the loading-1 face trace is -2 per unit load, giving `p_faces_0 = 2*(nu_u-nu)/[alpha*(1-2*nu)] = 0.7490636704`. The auxiliary states involving Omega, B, and beta are inconsistent. The loading-2 coupling factor `eta = alpha*(1-2*nu)/[2*(1-nu)] = 0.33375` is missing from the reported long-time results.

## 3. Method

**Elastic-limit assignments plus Cauchy-kernel quadrature; general_or_wrong.** The delivered program evaluates elementary expressions in `exact_results`, assigns constants in `coupled_results`, and separately tests a principal-value integration helper. It contains no coupled Biot solver, collocation matrix, time stepping, inversion, or short-time fit. The correct loading-1 volumes follow from the elastic opening; replacing the initial crack-tip solution by an invented elastic factor is invalid.

The write-up claims a Rice–Cleary reduction and Chebyshev boundary-integral solve (ANSWER.md:76–78, 118–123), but the delivered files do not demonstrate them. A full transient inversion is not necessary for the requested asymptotics. The key's original face-diffusion/reciprocity calculation also needs its stated tip correction: the face contribution alone is 1.354055, and the corrected full coefficient is 2.440220. Neither contribution is computed here.

## 4. Errors

Counts below refer to distinct errors, with repeated occurrences across prose, formulas, and code grouped together. Six erroneous requested outputs are counted individually as decisive. Further contradictory auxiliary equations, unsupported validation claims, and explanations accompanying correct quantities are counted as side errors. The inactive equations are not counted again as causes of the assigned numerical outputs.

### Decisive: 6

| ID | Evidence | Error and consequence |
|---|---|---|
| D1 | ANSWER.md:20–21, 36: “the incompressible-undrained limit raises K by 1/(1-nu)”; solve.py:83 | The initial loading-1 intensity is not `sqrt(pi*L)/(1-nu)`. This ignores the drained tip and gives 1.25 instead of 0.807593. It also fails the uncoupled limit `nu_u -> nu`, where the elastic intensity is 1, not `1/(1-nu)`. |
| D2 | ANSWER.md:40, 67, 117–123: “report 0.0; the true A is non-zero”; solve.py:118, 127 | The requested leading growth is missing, and zero is written as its coefficient. The corrected coefficient is 2.440220. The `formulas` entry merely defines A and repeats c; it supplies no expression for A. Flagging nonconvergence does not turn the placeholder into an answer. |
| D3 | ANSWER.md:125–129: “p = 1-sigma_n = 2”; results.json:17 | Pore pressure is incorrectly inferred from normal traction instead of undrained storage and stress trace. It must depend on alpha and both Poisson ratios. The supplied formula also says both `Omega = p` and `Omega = p+1`, and does not give the required general pressure expression. The result 2 should be 0.7490636704. |
| D4 | ANSWER.md:49–52: “the same as ... equivalent total traction”; solve.py:86 | A prescribed pore-pressure step with zero total face traction is incorrectly treated as an instantaneous mechanical traction. `K2_0` is 0, not -1.25. The finite negative intensity also conflicts with the submission's own explanation of zero initial displacement. |
| D5 | ANSWER.md:15–18, 52; solve.py:87 | Applying the material-independent dry mechanical intensity to loading 2 overlooks the maintained pore-pressure boundary condition and poroelastic back-stress. `K2_inf` is `-eta`, not -1. The general expression must contain alpha. |
| D6 | ANSWER.md:55, 71, 145–146; solve.py:88 | The long-time loading-2 volume is a closure, `-alpha*pi*(1-2*nu)*L^2/(2*G)`, not the positive loading-1 drained volume. The additional assumption `w >= 0` excludes the signed closure requested by this linear model. The reported value +2.5132741229 should be -0.8388052385. |

### Side: 11

| ID | Evidence | Error |
|---|---|---|
| S1 | ANSWER.md:23–28: “Loading 2 supplies no forcing to the Omega equation” | This contradicts the submission's definition `Omega = p-sigma_n`: loading 2 gives Omega = 1 and B = Omega+p = 2 on the faces. The claimed storage/force split is not established by the supplied equations. |
| S2 | ANSWER.md:125–127; solve.py:100: “Omega = -p, B = 0” | These asserted loading-1 initial states contradict its prescribed face data and the solver's own definitions, which give Omega = 1 and B = 1 at p = 0, sigma_n = -1. The opening and volume are nevertheless assigned the correct undrained elastic formula. |
| S3 | solve.py:9–11: “beta = sqrt(1-S) = (1-2nu)/(2(1-nu))” | The asserted equality is false: the two sides are 0.8660254038 and 0.375 for the stated parameters. |
| S4 | solve.py:13–16: “mass balance S dw/dt = (c/2) d2Omega/dx2” and “Biot solution” | This is not a valid reduction as stated. On loading-1 faces Omega is spatially constant, so this equation forces `w_t = 0`, contradicting the different correct initial and final volumes. It replaces bulk pore-fluid diffusion by a crack-only equation without a justified reduction. The following mechanical equations introduce an undefined auxiliary `y` and no demonstrated steady pressure-coupling dependence on alpha. |
| S5 | solve.py:13: “dw/dx = 0 at |x| = L” | Crack opening vanishes at the tips; the Sneddon opening's derivative diverges there. A zero tip slope contradicts the very opening profile used for the correct volumes. |
| S6 | solve.py:16: “Int_{-L}^{L} B dx = 0” | With the stated B definition and face data the integral is 2L for loading 1 and 4L for loading 2, not zero. No different B or separate constraint is supplied. |
| S7 | solve.py:18: “delta response (w=0 for t>0)”; ANSWER.md:28, 54, 130–131 | A pore-pressure step causes continuing diffusion and a nonzero long-time response, not an aperture identically zero for positive times. This also contradicts the submission's own nonzero `V2_inf`. Zero initial volume is correct, but neither this delta-response explanation nor the claim that it could be wrong only for a cubic-law conduit is justified. |
| S8 | ANSWER.md:76–78, 105–109, 118–123; solve.py:94–95 | The delivered code does not substantiate the claimed collocation solve, order-one residuals, nonconvergent inversion, or odd-Chebyshev test with worst difference `1.44e-15`. `coupled_results` contains only three assignments. This is a reproducibility/evidence gap; it does not establish that an unsubmitted historical experiment never happened. |
| S9 | ANSWER.md:89–90: “six-figure agreement rules [pole mishandling] out”; solve.py:45–51, 60–62 | For the constant-function check, `gfun(s)-gval` is identically zero; the remaining logarithm is precisely the expression used as the expected answer. Agreement cannot validate the quadrature or its pole treatment. A grader experiment replacing quadrature by a zero-returning function still passed Check 1. |
| S10 | ANSWER.md:101–103: “validates ... the drained value V ... used in (a) and (b)” | The dry-crack Hilbert identity verifies a kernel and an elastic profile, not the pressure-coupled loading-2 volume. Passing it does not justify applying the positive dry mechanical volume to loading 2. |
| S11 | ANSWER.md:138–144: “The problem does not state whether ... prescribed ... or ... cubic law” | The problem explicitly prescribes total face traction and pore pressure for both loadings. An unspecified internal conduit law is not missing information for this stated boundary-value problem. This weak ambiguity claim does not explain or support a different answer. |

### Presentation: 1

- **P1 — ANSWER.md:118–120:** calling the prescribed-loading poroelastic boundary-value problem “an inverse-scattering problem” is incorrect terminology; no inverse problem based on scattering observations is posed.

The repeated “exact” labels and Check 4's assertions are already covered by the erroneous results and reasoning above. Normalization and the stated K and V units are otherwise consistent.

## 5. Checks

**Four checks are reported; one demonstrably independent check could fail.** None tests the deciding poroelastic tip intensity, short-time coefficient, induced pressure, or loading-2 coupling. Individual sample points within a numbered check are not counted separately.

| Reported check | Tests the deciding quantity? | Could fail / evidence |
|---|---|---|
| 1. Constant-function PV kernel (ANSWER.md:80–90) | No | No independent quadrature test: the residual integrand is zero and both sides reuse the same logarithm. It still passes when the integrator is deliberately broken. Printed numbers reproduce. |
| 2. Dry-crack semicircle PV identity (ANSWER.md:92–103) | Only an elastic kernel/profile, not the deciding poroelastic results | Yes: numerical quadrature is compared with the independently known `-x`. It reproduces and fails when the integrator returns zero. Counted. |
| 3. Odd-Chebyshev weighted PV identity (ANSWER.md:105–109) | No; at most a spectral-operator identity | Potentially falsifiable if the integral is independently evaluated, but no implementation or full output is delivered. Comparing analytic `U` with its sine representation alone would be an identity. The claimed execution and accuracy cannot be credited as a demonstrated fail-capable test. |
| 4. “Elastic limits are analytic” (ANSWER.md:111–113) | It asserts the disputed intensities | No: there is no residual, independently computed comparison, acceptance condition, or limit experiment. Analytic assignment and independence from discretization do not establish correctness. |

The grader independently confirmed that the odd-Chebyshev mathematical identity itself is valid; this does not reproduce the solver's claimed run. `checks_could_fail` counts the independent tests demonstrated by the delivered material, hence 1. Zero of those checks validates the decisive poroelastic quantities. The two principle checks in AUTO_CHECKS.json are grader checks, not additional solver-reported checks.

## 6. Claims

**Partly supported overall.** AUTO_CHECKS.json records exit 0 and all nine quantities printed. A fresh run of `python3 solve.py` on a copy in `scratch/reproduction` also exited 0; the rewritten results.json exactly equals the submitted JSON, including formulas. All nine numerical claims and the pasted outputs of Checks 1 and 2 reproduce. The loading-1 volume ratio is 4/3 as required, and `K1_inf = 1` passes the Griffith check. Reproducibility does not make the six failed quantities correct.

The claimed coupled numerical implementation and Check 3's specific precision/output are absent. Check 4 has no executable verification. Thus numerical-output agreement is complete, but support for the full set of claimed computations and checks is only partial. No numerical reproduction mismatch was found. Detailed grader evidence is saved in `scratch/numeric_audit.json`, `scratch/check_strength.txt`, and `scratch/additional_audit.txt`.

## 7. Questions

**Four numbered questions** appear in ANSWER.md:136–150, plus a separately flagged failure to determine A.

1. The crack-model question is unnecessary: prescribed face tractions and pressure already define this problem. A cubic-law conduit would be a different, incompletely specified model.
2. The V-definition question repeats the stated opening integral and unit thickness; its added nonnegative-opening assumption is inappropriate for loading-2 closure.
3. The K-normalization question uses the correct classical convention and agrees with the reporting units, but identifies no missing information.
4. Distinguishing pre-diffusion induced pressure from the prescribed zero face pressure is sensible. Equating that pressure to Omega is not justified and contradicts the solver's definition, as D3 records.

Admitting uncertainty about A is candid, but missing convergence of an unsubmitted calculation is not missing problem data.

## 8. Source

**Authors are named.** solve.py:6–7 states: “Rice & Cleary 1976; Detournay & Cheng 1993”. ANSWER.md:142 also names “Rice–Cleary/Detournay”. No paper title is given. The dated Detournay–Cheng reference does not identify the key's 1991 crack paper; without a title, it cannot be determined whether the solver intended that paper or a different work. Author-naming credit is warranted; no unsupported assertion about the identity of the 1993 source is needed.

## 9. Challenges

**One weak ambiguity challenge** is present (ANSWER.md:138–144): the solver claims the crack model is unspecified and a fluid-filled crack governed by the cubic law would change the answer. The argument only contrasts prescribed face conditions with an alternative conduit model and mentions `h^3`; it supplies no alternative initial/boundary conditions, parameters, computed result, or limit for another person to reproduce. It is **not checkable as presented**, and the given problem already specifies the face conditions. It does not qualify for the well-argued, checkable-challenge exception; S11 records the erroneous ambiguity assertion.

There is no challenge to the corrected initial tip intensity, the tip contribution to small-time growth, or a published numerical result. No answer is contested. This audit finds no evidence requiring an additional correction to KEY.md; it applies the explicit P374 corrections rather than the obsolete values retained for historical context.

```json
{"answers": {"a": "partly", "b": "wrong", "c": "partly"},
 "supporting_state": "wrong",
 "method": "Elastic-limit assignments plus Cauchy-kernel quadrature", "method_class": "general_or_wrong",
 "errors": {"decisive": 6, "side": 11, "presentation": 1},
 "checks_reported": 4, "checks_could_fail": 1,
 "claims_supported": "partly",
 "questions_listed": 4,
 "source_named": "Rice & Cleary 1976; Detournay & Cheng 1993",
 "challenges": [{"claim": "Crack model is unspecified; a cubic-law conduit would change the answer", "checkable": false}]}
```
