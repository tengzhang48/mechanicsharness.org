# Grading profile

## 1. Answers

The corrected key at KEY.md:11–34 controls grading; its later historical paper formulas do not override that correction. All nine numerical results pass the stated bands. The contested labels below record well-supported interpretation challenges, not numerical errors.

- **(a) Contested; right against the corrected key, right under the solver's explicit limit convention.** The true-tip K1_0 is 0.807593003, K1_inf is 1, the volume endpoints are π(1−ν_u)L²/G and π(1−ν)L²/G, and the complete short-time coefficient is 2.440219579; the reverse-order limit K1_0=1 is also explicitly supplied (ANSWER.md:36–65, 321–326).
- **(b) Right.** K2_0=0, K2_inf=−α(1−2ν)√(πL)/[2(1−ν)], and V2_inf=−απ(1−2ν)L²/(2G) have the correct magnitudes and closure signs, with a coherent transient long-time argument (ANSWER.md:67–76, 288–315).
- **(c) Contested; all nine values right against the key.** Every requested number and its general formula is delivered; the solver also gives the alternative instantaneous outer K1_0=1 and literal face-pressure trace 0, so its answer remains correct under the stated alternative readings (ANSWER.md:19–34, 317–338; results.json).

## 2. Supporting state

**Right.** The initially undrained ellipse, its face stress trace s=−2, induced outer pressure 2(ν_u−ν)/[α(1−2ν)]=0.7490636704, drained/undrained volume ratio 4/3, and Griffith long-time intensity are correct. The actual prescribed pressure trace is correctly distinguished from this outer pressure (ANSWER.md:105–128, 327–330).

The solver locally defines η=(ν_u−ν)/(1−ν_u)=1/3; the key uses the same letter for the different back-stress factor 0.33375. The solver correctly evaluates that back-stress separately. Its factor entries n11=1.07679067033231 and n21=0.416959452647314 give n11²+n21²−(1+η)=−4.67×10⁻¹⁰. The face contribution 1.354055001 and tip contribution 1.086164579 sum to the reported 2.440219579. The harmonic variable h, fluid-content coefficient D_f, and half-space mixed boundary conditions are consistent with the supplied equations (ANSWER.md:107–161).

## 3. Method

**Right method:** bulk drained/undrained elasticity, Fourier–Laplace surface reduction, canonical matrix Wiener–Hopf factorization of the local tip problem, and undrained elastic reciprocity for the volume change. This resolves the drained inner tip and includes the same-order tip contribution omitted by face diffusion alone. It addresses precisely the two quantities corrected in the key; it is not an unnecessary full transient computation.

The delivered primary code actually factors the matrix on successively refined mapped frequency grids and extrapolates its entries; the two deciding numbers are not hard-coded (factor_check.py:7–49; solve.py:17–32). The finite-crack Fourier collocation/CG solve is a separate numerical audit, using the coupled boundary symbol but no factored amplitudes as inputs (finite_check.py:14–87).

## 4. Errors

**No submission errors found:** decisive 0, side 0, presentation 0. The limit qualifications are justified, the negative loading-2 values follow the stated linear model, and the numerical uncertainty and unsuccessful refinements are disclosed. The finite-crack check does not certify all displayed digits, and the solver explicitly says so (ANSWER.md:353–360; finite_check.md:95–136). The differing definition of η is explicit and consistent.

**Key transcription issue, not a submission error:** KEY.md:61 writes the opening–K prefactor as √(πG)/[√8(1−ν)]. Dimensional consistency and substitution of w=2(1−ν)√(L²−x²)/G require G√π/[√8(1−ν)] instead. The printed prefactor would give K=√(πL/G) for unit traction. Against that literal prefactor, the submission's general G dependence disagrees; against the dimensionally corrected relation and the key's explicit normalized targets, its intensities are right. This ancillary transcription error does not override those targets. The key's obsolete 0.75 and 1.354055 results are expressly superseded by its correction.

## 5. Checks

Counting convention: four numerical audit sections in CHECKS.md:7–161; the four distinct checks in its section 5; and the eight finite-crack campaigns in finite_check.md:159–421. This gives **16 reported check entries**, including one explicitly incomplete attempt. Repeated output in the finite-check summary and individual grid points are not counted again. **13 completed entries could reject a wrong mathematical result.** Reproduction and schema checks can reject delivery defects but do not test the deciding physics; the timed-out campaign does not supply its intended completed comparison.

| Entry and reference | Deciding quantity tested | Could reject a wrong mathematical result? |
| --- | --- | --- |
| P1: initial refinement, CHECKS.md:7–28 | K1_0 and hydraulic factor n21 through refinement; factor norm, reality, and projection residual | Yes, for factorization/resolution errors; a residual alone cannot validate the governing symbol. |
| P2: extended refinement, CHECKS.md:30–46 | Same amplitudes on the additional 1,048,576-point grid | Yes; repeats the numerical method at finer resolution, with no independent model validation. |
| P3: frequency-map scale comparison, CHECKS.md:48–75 | Extrapolated K1_0 and A across three map scales | Yes; can expose frequency clustering and extrapolation bias. |
| P4: final extrapolation, CHECKS.md:77–161 | Raw factors, linear/quadratic extrapolated K1_0 and A, and norm identity | Yes; checks numerical stability and the coupled amplitudes, with shared-model limitations. |
| P5: exact reproduction, CHECKS.md:163–179; check_results.py:12–23 | Stored, printed, and regenerated results | No physical test; can fail if those deliveries disagree. |
| P6: schema, check_results.py:24–25 | Number/string types of requested fields | No physical test; can fail for missing fields or wrong types. |
| P7: energy-derived coefficient, check_results.py:29–39 | A reconstructed from K versus reported A | Yes; a hydraulic-amplitude, sign, or coefficient error can violate the relation. This is an identity check using the same primary factorization. |
| P8: weak coupling, check_results.py:41–45 | (1−K)/η approaches 3/4 at η=10⁻⁴ | Yes; tests the tip method's limiting behavior, not all finite-coupling answers. |
| F1: initial finite-crack run, finite_check.md:159–198 | A from volumes across grids and Laplace parameters; same-grid versus exact-reference subtraction | Yes; directly tests volume growth and reveals unresolved layers/reference bias. |
| F2: larger grids and higher parameters, finite_check.md:200–233 | Same A on finer grids, with high-parameter turnover | Yes; detects tip underresolution; coarse turnover is not accepted as the physical limit. |
| F3: box doubling, finite_check.md:235–257 | A at matched spacing in boxes 32 and 64 | Yes; tests periodic-box dependence, with remaining grid/parameter errors. |
| F4: million-point refinement, finite_check.md:259–287 | A on two finer grids through parameter 10⁸ | Yes; directly exposes remaining discretization dependence. |
| F5: box halving, finite_check.md:289–310 | A at matched spacing in boxes 16 and 32 | Yes; can expose a biased infinite-domain estimate. |
| F6: timed-out finest run, finite_check.md:317–348 | Intended finer-grid/high-parameter comparison | Incomplete; partial volumes/residuals exist, but no completed intended limit check is credited. |
| F7: completed box-64 extrapolation, finite_check.md:350–384 | Direct volume-based A, grid extrapolation, and alternative large-parameter fits | Yes; estimates 2.438239292 and 2.440642227 could disagree with the primary result. The three-point quadratic residual is not an accuracy bound. |
| F8: matched-spacing box-32 run, finite_check.md:386–421 | Same volume coefficient and fits, compared across boxes | Yes; tests box sensitivity of the limiting estimate. Assumed box-error powers remain empirical. |

These are falsifiable comparisons, not thirteen independent derivations. Factor checks share a kernel; finite-crack checks share a surface reduction. The finite-crack solver integrates the opening directly and therefore can detect omission of the tip volume term. It does not extract K, simulate loading 2 at long times, or perform time stepping; those checks are not claimed.

## 6. Claims

**Supported: yes.** AUTO_CHECKS.json reports successful reproduction of all nine printed claims and passes all corrected-key and principle checks. On an untouched copy in scratch/submission/, rerunning check_results.py exited 0 and reproduced its reported output, including exact JSON equality, energy-coefficient difference 1.458×10⁻⁹, and weak-coupling deficit 0.749865248. The regenerated results also exactly match the original delivered JSON. Logs are in scratch/check_results_rerun.log and scratch/review_identities.json.

The reported box-64 finite-crack command also exited 0 on the copy. Its volumes, coefficients, residuals, extrapolations, and fit results reproduce the transcript; elapsed run times naturally differ (scratch/finite_check_rerun.log). No numerical claim mismatch was found. Earlier-version factorization and finite-check commands are identified as historical in the submission, so the current script's additional output is not a concealed discrepancy. The transcript's timed-out run is expressly not claimed as completed.

## 7. Questions

**Six listed entries** under Questions (ANSWER.md:317–338): unstated shear traction; K1_0 order of limits; induced outer pressure versus prescribed face trace; the two-dimensional loading-2 long-time limit; absence of contact constraints; and the ν_u=1/2 endpoint/general factor-defined formulas. The last is a mathematical qualification rather than a missing datum. The questions and conventions are sensible. Zero shear is the usual symmetric mode-I assumption; signed closure follows the supplied linear equations. No answer is withheld for missing information.

## 8. Source

**None.** The submission names no source paper or author for this crack problem. “Canonical Wiener-Hopf factorization” (ANSWER.md:97) names a method, not an attribution to the source paper. The source named by the key is not credited to the solver.

## 9. Challenges

1. **Face-only short-time growth is incomplete.** Quote: “a planar diffusion layer alone misses an equally large tip contribution” (ANSWER.md:101–103). Argument: the singular auxiliary stress makes tip drainage contribute at the same √t order; the reciprocal/Green-identity calculation gives A_faces+A_tips (ANSWER.md:231–274). **Checkable:** factor the supplied matrix and independently solve the finite-crack mixed problem; the rerun supports A≈2.44. Against the corrected key, (a) and (c) are right; against the challenged face-only approximation, the solver's complete coefficient is supported. The old 1.354055 approximation is incomplete, not a solver error.
2. **K1_0 depends on the order of singular limits.** Quote: “K1_0 is ambiguous without an order of limits” (ANSWER.md:321). Argument: the outer initial response is undrained, whereas the inner tip at each positive time is drained; the coupled tip factor yields 0.807593003, while taking time to zero first gives 1 (ANSWER.md:12–15, 194–221, 321–326). **Checkable:** extract the tip intensity at positive times before extrapolating, or instead take the undrained outer solution first. The corrected key explicitly recognizes both readings. Against its selected tip convention, (a)/(c) are right; under the opposite convention, the solver explicitly supplies the right alternative. These answers are marked contested.
3. **The requested induced face-adjacent pressure differs from the boundary trace.** Quote: “The actual loading-1 face pressure trace is exactly zero” (ANSWER.md:327). Argument: the prescribed Dirichlet value is 0, while the initially undrained solid outside the vanishing diffusion layer has pressure 2d/a (ANSWER.md:61–65, 247–251, 327–330). **Checkable:** evaluate the face boundary value and the fixed-distance t→0 outer limit separately. Against the key's induced-pressure convention, 0.7490636704 is right; against a literal face-trace reading, the solver's explicitly stated alternative 0 is right. This contributes to the contested designation for (c).
4. **Loading 2 has no decaying classical two-dimensional steady pressure field.** Quote: “there is no decaying harmonic pressure field with p=1 on a finite slit in an infinite two-dimensional domain” (ANSWER.md:295–297). Argument: an expanding pressure halo gives a nonuniform spatial limit and selects the local back-stress −b/2 (ANSWER.md:295–315). **Checkable:** solve in increasing outer domains or take the small-Laplace-parameter limit with the far-field matching retained. This challenges a globally steady interpretation, while the problem asks for the transient long-time limit. The submission and key give the same requested limits, so (b) is right under both and has no numerical dispute.

The shear/contact/endpoint remarks state assumptions or possible model extensions; they do not assert that the stated linear Biot problem is wrong. No source paper or published result is explicitly accused by name. None of the four supported, checkable challenges above is counted as an error.

```json
{"answers": {"a": "contested", "b": "right", "c": "contested"},
 "supporting_state": "right",
 "method": "Bulk elasticity, Fourier-Laplace/Wiener-Hopf tip factorization, reciprocity",
 "method_class": "right",
 "errors": {"decisive": 0, "side": 0, "presentation": 0},
 "checks_reported": 16, "checks_could_fail": 13,
 "claims_supported": "yes",
 "questions_listed": 6,
 "source_named": "none",
 "challenges": [{"claim": "Face-only short-time volume growth omits a same-order tip contribution", "checkable": true}, {"claim": "Instantaneous stress intensity depends on the order of tip and time limits", "checkable": true}, {"claim": "Induced outer pressure differs from the prescribed face trace", "checkable": true}, {"claim": "Loading 2 has no decaying classical two-dimensional steady pressure field", "checkable": true}]}
```
