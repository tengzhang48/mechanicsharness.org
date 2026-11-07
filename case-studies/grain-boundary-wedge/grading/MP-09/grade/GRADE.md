# Grading profile

## 1. Answers

- **(a) Right.** The computed eigenproblem, initial transient, and first-mode exponential tail support the requested evolution; lambda1 = 11.7426682103 and lambda2 = 264.6342256450 meet the key's 1% and 2% bands. The code shows the early overshoot that the answer describes; its unqualified word “monotonically” needs a late-time qualification.
- **(b) Right.** The limit is a traction-free, crack-like plated wedge, and the computed K_inf_norm = 1.1215223627 lies in 1.09–1.15.
- **(c) Right.** relaxation_time_t0 = 0.0851595210 equals 1/lambda1 and lies in 0.05–0.15.

These judgments credit the actual computations and the correctly described transient and tail. They do not endorse every explanatory statement or claimed precision.

## 2. Supporting state

**Right for the quantities requested by the problem.** The slowest-mode profile is normalized to f(1) = 1, and f(0.5) = 0.8019565032 agrees with the key's 0.8020 within 2%. The production basis satisfies f(0) = 0 and f'(1) = f''(1) = 0. The optional amplitude and statements about endpoint asymptotics have errors listed below; they are not requested supporting quantities.

## 3. Method

**Right: singular-integral spectral collocation with SVD reduction.** Combining diffusion and elasticity gives sigma_tau = A sigma and A f = -lambda f with the correct sign, scaling, and physical boundary conditions. The s = sqrt(1-x) basis resolves fractional tip powers; analytic principal-value moments and numerical quadrature for the regular image terms form the collocation pencil. Truncated SVD addresses its numerical singularity. This is the singular-integral numerical method the key calls for, with a different basis. The finite-difference discretizations provide approximate cross-checks. The steady crack integral equation and the K/COD normalization are correctly implemented.

## 4. Errors

**Decisive: 0.** No identified error changes the requested rates, profile value, crack limit, or relaxation scale. Repeated statements of the same error are grouped below.

**Side: 12.**

1. **Wrong leading tip behavior.** `submission/ANSWER.md:28–29` says “eigenfunctions behave as (1−x)^(7/2)” and “f''' ~ (1−x)^(1/2)”; `submission/sblib.py:3–10` repeats this. The normalized mode tends to 1, with leading correction proportional to (1-x)^(5/2), and f''' proportional to (1-x)^(-1/2). Its computed s^5 coefficient is -1.2251213; sqrt(1-x) f''' tends to 2.29710. The code includes s^5 and therefore captures the correct behavior despite this explanation.
2. **Wrong kernel statement.** “ker(A) = span{1}” (`submission/ANSWER.md:26`) is false for the unconstrained differential-integral expression: A annihilates 1, x, and x^2. The three imposed boundary conditions exclude all nonzero members of that polynomial subspace, not just the constant.
3. **Incorrect translation of the tip conditions.** `submission/sblib.py:8` gives f_s(0) = f_ss(0) = f_sss(0) = 0 and says this “kills the ker(A) mode s^4.” But s^4 satisfies all three of those conditions while f_xx(1) = 2. The actual basis correctly excludes s^4; the stated equivalence and justification are wrong.
4. **Overstated precision and stability.** “accurate to ~5e-7” for lambda1, lambda2, and f(0.5), and “lambda1 then varies < 5e-7” (`submission/ANSWER.md:45–46,97,166–167`; `submission/sblib.py:110`) are not supported as absolute bounds. The delivered N=60 and N=110 lambda2 values differ by 4.67e-5. A scan over the stated resolutions and SVD tolerances gives a lambda1 spread about 7.53e-7 and lambda2 spread about 9.83e-5. Convergence spread also does not establish an error bound against the continuous solution. All values remain well inside grading bands.
5. **False Richardson bracket.** `submission/solve.py:11–12,235–236,275–276` and `submission/results.json:14` claim that the Richardson limits bracket the spectral rate. The supplied log gives [11.74267930, 11.74268899], with spectral 11.74266821 below both (`submission/run.log:19–20`). The raw finite-grid rates do bracket it; the extrapolated values do not. The answer partially acknowledges the extrapolation error, but the claimed bracket remains false.
6. **Misread FD operator validation.** “the interior error is O(dx)” (`submission/ANSWER.md:127–130`) does not describe the reported four-cell-trimmed maximum. It stays near 0.191 for FD-A and 0.237 for FD-B as N increases. The maximizing locations move toward the mouth; excluding four cells does not leave a fixed physical interior. A separate fixed interval [0.05,0.95] does converge, so this is an interpretation error rather than evidence that the primary solution is wrong.
7. **False conditioning claim for the crack solve.** “well-conditioned (c = O(1))” (`submission/ANSWER.md:107`; `submission/sblib.py:163`) is contradicted by the delivered monomial matrix. Its 2-norm condition number is about 4.3e8 at M=10 and reaches approximately 1e17–1e18 at M=18–24; coefficients reach about 656. Stable reported K values do not make the matrix well-conditioned.
8. **Overstated independence of K/COD checks.** The two K routes “only agree if the solve is correct” (`submission/ANSWER.md:159–162`) is false. For any density with tip coefficient c0, the outside-stress limit gives -pi sqrt(2)c0, whether or not the traction equation was solved. An intentionally wrong c0=-0.1 gives matching K limits approaching 0.444288 while its traction integral is only about 0.35–0.46 instead of 1. The near-tip COD–K relation likewise follows from the same c0. These are useful convention/implementation checks, not independent validation of the crack solution; comparison with the handbook value is a separate substantive check.
9. **Biased optional amplitude labeled as an asymptotic coefficient.** `submission/results.json:7` gives c1=1.54733364, and `submission/ANSWER.md:171–177` calls the quoted estimate the “honest late-asymptotic amplitude.” `submission/solve.py:175–189` instead back-extrapolates a finite-time FD/Euler solution using the different spectral rate. Biorthogonal projection on the same FD operator gives c1=1.5307202 at N=512 and an extrapolation near 1.53105. The growing estimate also includes Euler rate bias and snapshot timing effects, beyond the stated FD eigenrate offset. The coarse “approximately 1.5” description remains reasonable.
10. **Initial modal projection incorrectly ruled out.** “a purely modal projection of the initial condition is not well-defined at t = 0” because the initial trace violates the mouth condition (`submission/ANSWER.md:176–177`) does not follow. A pointwise boundary incompatibility permits a weak/mild initial-value solution and a modal expansion in the interior; it does not prevent modal coefficients. Left/right eigenvector projection of uniform initial data is well defined and converges in the delivered FD discretizations.
11. **Wrong expected mouth limit in an auxiliary test.** “R_m(s->1) -> 0 (kernel vanishes at x=0)” (`submission/test_slib.py:23–26`) incorrectly interchanges a nonuniform limit and integration. In particular, R_1(s) = T_0(1-s^2)/2 tends to -2. The rerun prints a maximum near 2, directly contradicting the expected zero. The main solver does not impose this incorrect expectation.
12. **Obsolete, unreliable auxiliary eigenvalue check.** `submission/test_slib.py:43–57` still uses the unreduced pencil instead of the production SVD fix. At N=44 its purported lambda1 is 388.23938 and f(0.5)=1.01075. It is unsuitable as a successful convergence check; completion of the script does not validate those numbers. The corrected production route supplies the graded results.

**Presentation: 7.**

1. “The traction relaxes monotonically” (`submission/ANSWER.md:42`) needs “at late times.” Its own transient and output show traction rising above the initial value. The correctly described overshoot and computed late tail resolve the intended evolution, so this wording does not change the answer grade.
2. “sigma(0,0+) = 1” (`submission/ANSWER.md:34`) mislabels the incompatible initial trace: for every positive time the mouth boundary value is zero. The initial uniform datum is 1 in the interior.
3. “c1 approximately 1.53 sigma0” (`submission/ANSWER.md:171`) attaches stress units to the dimensionless coefficient already multiplying sigma0 in line 39.
4. “t_relax = 1/lambda1 = 0.085160 t0” (`submission/ANSWER.md:71`) omits t0 in the first equality. The dimensional expression is t_relax = t0/lambda1; the JSON convention is correct.
5. “The kernel action is exact” (`submission/ANSWER.md:85–92`) overstates the implementation: the singular moments are analytic, but the regular image contribution uses finite Gauss–Legendre quadrature (`submission/sblib.py:39–55`).
6. The process log writes “BCs sigma(0)=0, sigma(1)=sigma(1)=0” (`submission/PROCESS.md:1`), omitting the derivative marks. The answer and production code use the correct conditions.
7. “positive rates” as a cause of non-decay (`submission/ANSWER.md:148–149`) confuses growth eigenvalues of A with the defined decay rates lambda. Positive lambda gives decay; positive eigenvalues of A give growth.

The earlier K post-processing mistake described as fixed in the submission is not counted as a remaining error: inspection and reproduction confirm the correction.

## 5. Checks

For transparent counting, the five bundled paragraphs in the answer are separated into **13 individual comparisons**, and the two auxiliary scripts contain **9 labeled test blocks**: **22 reported checks**. “Could fail” means the displayed numerical comparison can disagree; it does not imply a passing, independent, or sound correctness test. **21 could fail** under that definition. None of the scripts enforces these comparisons with pass/fail assertions.

| Answer check | Tests the deciding quantity? | Could fail? / finding |
|---|---|---|
| FD operator versus analytic Af | Indirectly: discretized governing operator | Yes; the reported trimmed-error convergence claim fails, although a fixed physical interior converges. |
| Spectral resolution convergence | Yes: both rates and requested profile | Yes; supports the grading bands, not the claimed absolute precision. |
| SVD tolerance sensitivity | Yes: rate stability | Yes; not printed by the reproduction command; a separate scan exceeds the stated bound. |
| Off-grid eigen-equation residual | Yes for the first mode | Yes; about 1e-6 relative at four noncollocated points. No corresponding second-mode residual is provided. |
| FD rate convergence / Richardson agreement and bracketing | Yes: both rates, chiefly lambda1 | Yes; raw-grid bracketing succeeds, extrapolated bracketing fails. |
| FD first-mode half-profile comparison | Yes: mode1_shape_half | Yes; the two discretizations approach and bracket the spectral value. |
| Time-domain late decay fit | Yes: relaxation rate | Yes; 11.712747 is close to the FD/Euler rate 11.711764, not the exact spectral rate. It checks evolution of the FD operator. |
| Time-domain decay toward zero | Yes: limiting traction | Yes; final maximum 0.003853 is still decaying. A finite endpoint alone is not proof of the infinite-time limit. |
| Crack M-convergence | Yes: K and mouth opening | Yes; K_tip stabilizes near 1.121522 despite ill-conditioning. |
| Tip versus outside-stress K | K extraction/conventions, not independent crack-equation correctness | Yes for an extraction implementation error; the reported earlier basis bug illustrates this. A wrong density can pass. |
| Computed K versus classical 1.1215 | Yes: part (b) normalization and value | Yes; the computed value agrees. The reference number itself is merely printed as a constant. |
| Computed mouth opening versus 5.8172 | Independent reference for an optional state quantity | Yes; computed 5.81719168 agrees. |
| Near-tip COD–K relation | No independent test of the solved density | No as delivered: no separate calculation/output is supplied, and the leading relation uses the same c0 on both sides. |

| Auxiliary check | Tests the deciding quantity? | Could fail? / finding |
|---|---|---|
| `test_lib.py` 1: infinite-domain kernel moment -4 | Indirectly: kernel/sign normalization | Yes; reproduces -4. |
| `test_lib.py` 2: monomial moments versus PV quadrature | Indirectly: operator moments | Yes; agrees to about 1e-14. |
| `test_lib.py` 3: T0 versus antiderivative | Indirectly: operator moments | Yes; agrees to roundoff. |
| `test_lib.py` 4: Chebyshev action versus analytic Af | Indirectly: spectral operator assembly | Yes; agrees to about 1e-13. |
| `test_lib.py` 5: integrated-basis derivatives and BCs | Indirectly: basis constraints | Yes; finite differences are consistent with the expected derivative/BC limits. |
| `test_slib.py` S1: s moments versus quadrature | Indirectly: regular image moments | Yes; largely agrees. Both routes share the same analytic Cauchy moment, limiting independence. |
| `test_slib.py` S2: zero mouth limit | Wrong expected limit, not a sound correctness test | Yes, and visibly fails: output is about 2, not 0. Counted as falsifiable, not as successful validation. |
| `test_slib.py` S3: s operator versus analytic Af | Indirectly: production moment machinery | Yes; agrees, with near-mouth numerical error around 4e-9. |
| `test_slib.py` S4: unreduced eigenvalue convergence | Intended to test rates/profile, but does not validate the corrected solver | Yes, and fails badly at larger N. |

## 6. Claims and reproduction

**Partly supported overall; all required numerical fields supported.** `AUTO_CHECKS.json` reports exit status 0, all 9 numerical claims printed, exact-class lambda1, and passing lambda2, profile, K, time, ordering, and reciprocal-time checks. A new `python3 solve.py` run on a copy in `scratch/reproduction/` exits 0 and reproduces every required numeric field exactly. The optional c1 differs by only 7.8e-10 after limiting BLAS threads, which is immaterial roundoff rather than a substantive reproduction mismatch.

The logs reproduce the advertised FD errors, spectral residuals, time decay, and K routes, including the failed Richardson bracket. They do not support the claimed universal precision, well-conditioning, correct auxiliary endpoint limit, or successful unreduced eigenvalue convergence. The near-tip COD comparison “coefficients 6.3435 vs 6.3443” (`submission/ANSWER.md:161–162`) is not produced by any delivered script. A SVD-tolerance scan is also claimed but not part of the advertised reproduction output. Printing reference constants and finishing a script do not establish that its comparisons passed.

Evidence is retained in `scratch/reproduction/grader-rerun.log`, `manual-review.log`, `amplitude-review.log`, `test-lib-rerun.log`, and `test-slib-rerun.log`. Original submission files were not changed.

## 7. Questions

**0 questions or missing-information requests.** “None blocking” (`submission/ANSWER.md:183`) is followed by sensible convention notes about 2u, COD, the K normalization, and plane-strain E*. The given model supplies the information needed for the requested answers.

## 8. Source

The only named author is **Tada**, quoted as “vs classical 1.1215 (Tada)” (`submission/ANSWER.md:156`), with the mouth-opening reference also attributed to Tada in the run output. This is a crack-handbook attribution, not identification of the grain-boundary-diffusion source. The submission does not name Gao, Zhang, Nix, Thompson, Arzt, or the paper identified by the key.

## 9. Challenges

**None.** The submission discusses numerical conditioning, a repaired implementation bug, and initial-condition incompatibility; it does not argue that the problem, the key's model, or a published result is wrong. There is no contested answer. This review finds no evidence of a key error affecting the stated grading bands.

```json
{"answers": {"a": "right", "b": "right", "c": "right"},
 "supporting_state": "right",
 "method": "Singular-integral spectral collocation with truncated-SVD eigenanalysis; steady crack collocation",
 "method_class": "right",
 "errors": {"decisive": 0, "side": 12, "presentation": 7},
 "checks_reported": 22, "checks_could_fail": 21,
 "claims_supported": "partly",
 "questions_listed": 0,
 "source_named": "vs classical 1.1215 (Tada)",
 "challenges": []}
```
