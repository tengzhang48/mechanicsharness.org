# Submission grade

## 1. Answers

Graded against the **corrected P374 key**, using its 1% bands for nonzero numbers and absolute 0.01 band for `K2_0`.

- **(a) Partly right.** Initial/final volumes and final intensity are right; the provisional drained-tip intensity is also within tolerance. The leading short-time volume growth omits the tip contribution and is wrong. The general initial-intensity expression is an implicit mixed-operator definition, with its numerical uncertainty openly acknowledged.
- **(b) Right.** Both intensities and the final signed volume match the key, including the negative signs for closure.
- **(c) Partly right.** Eight of nine numerical entries pass; `V1_short_coeff` is 44.51% below the corrected key.

| Quantity | Submission | Corrected key | Judgment |
|---|---:|---:|---|
| `K1_0` | 0.8069769387 | 0.807593 | Right; relative error 0.0763% |
| `K1_inf` | 1 | 1 | Right |
| `V1_0` | 1.8849555922 | 1.8849555922 | Right |
| `V1_inf` | 2.5132741229 | 2.5132741229 | Right |
| `V1_short_coeff` | 1.3540550005 | 2.440220 | Wrong; missing approximately 1.086165 |
| `p_faces_0` | 0.7490636704 | 0.7490636704 | Right in the requested initial outer limit |
| `K2_0` | 0 | 0 | Right |
| `K2_inf` | -0.33375 | -0.33375 | Right |
| `V2_inf` | -0.8388052385 | -0.8388052385 | Right |

The uncorrected paper/key values would reverse the judgments on the two sensitive entries: 0.806977 would fail against 0.75, whereas 1.354055 would pass. Those superseded values are not used for this grade.

## 2. Supporting state

**Right.** The explicitly requested supporting quantity, the initial induced outer face-adjacent pressure, is correct. The solver appropriately distinguishes this from the exact boundary trace `p=0` at positive time (`submission/ANSWER.md:23`, `:43`, `:57`). The full opening, its integral, and the drained/undrained modulus choices are correct. Its definitions give `d=0.2`, `b=0.3745318352`, `eta=0.33375`, `r=0.75`, and `beta=0.5`. The auxiliary flux-label error below does not change this requested state.

## 3. Method

**Mixed Fourier–Laplace half-space operator, periodic collocation, elastic calibration, and empirical tip extrapolation; valid alternative.** The surface formulation enforces drained faces and symmetry outside the slit and resolves the drained tip separately from the undrained outer field. This is appropriate for the corrected initial-intensity question. The ellipse calculations suffice for the volume endpoints; the small-Laplace-frequency limit gives the correct loading-2 results.

The full periodic crack calculation is more work than a local tip factorization, but is a valid route, rather than an intrinsically wrong method. The extrapolation remains provisional: `solve.py` fits stored values, uses an estimated mesh correction at `ell=0.025`, and assumes a linear dependence on `sqrt(ell)`. Those limitations are disclosed. The short-time volume derivation incorrectly integrates a nonuniform interior expansion without the necessary tip contribution.

## 4. Errors

**One decisive error, two side errors, no presentation errors.** Repeated appearances of the same defective short-time formula are counted as one error.

1. **Decisive — incomplete short-time volume asymptotic.** `submission/ANSWER.md:40–42`, `:69–76`, `submission/results.json:6`, `:16`, and `submission/solve.py:19` report the face-only coefficient `16*d*(1-nu_u)/[(1-nu)*sqrt(pi)] = 1.3540550005`. The reasoning explicitly says, “The displacement expansion is not uniform at the endpoints, but its leading integral is the stated correction” (`ANSWER.md:76`). That last inference is false: the endpoint contribution also enters at order `sqrt(t)`. The corrected coefficient is 2.440220, including approximately 1.086165 from the tips. This changes (a) and (c), even though it agrees with the superseded paper formula.
2. **Side — flux derivative has the wrong stated sign.** `submission/ANSWER.md:80` defines “q=p_y,” but the second column of its matrix (`:83–84`) and `submission/tip_check.py:17–18` use `q=-p_y`. From the stated surface relation at `ANSWER.md:64`, if `q=p_y`, then `p0 = 2*b*Du*k*(m-k)*v/m - q/m`; the written matrix instead has `+q/m`. Its traction-column sign is reversed as well. Changing the auxiliary variable to `-p_y` reconciles both columns and preserves every displacement, volume, and intensity, so this is a mathematical convention error without an answer change.
3. **Side — misleading volume-slope validation.** `submission/ANSWER.md:157` says, “The volume divided differences are consistent with the predicted Laplace short-time slope 1.2.” Its own fixed-`n=2048` data give `(1.9794399542-1.9422132554)/0.025 = 1.48907` and `(1.9422132554-1.9213397890)/0.0125 = 1.66988`, with an increasing discrepancy as the diffusion length decreases. These are not validating evidence for 1.2. The corrected key corresponds to a Laplace slope `2.440220*sqrt(pi)/2 = 2.16259`. The incorrect assertion is distinct from the missing-tip derivation; the solver says the comparison was not used for precision validation or its final fit.

The doubled volumes captured at `submission/ANSWER.md:103–107` are expressly **rejected** at `:112` and corrected in the delivered answer, JSON, and code. They are recovered mistakes, not additional surviving errors. Likewise, the candid uncertainty on the initial intensity is a limitation, not an erroneous claim of exactness.

## 5. Checks

Count **16 reported check types**, grouping repeated frequency/grid samples of one test together. A check counts as capable of failing when it can reject an algebraic, normalization, numerical, or asymptotic error in the quantity it actually compares. This does not imply it rejects every wrong physical solution. Delivery-only tests are identified separately and are not counted as solution checks. **Nine check types could fail on a solution or operator error.**

| # | Reported check and location in `submission/ANSWER.md` | Tests the deciding tip intensity or growth coefficient? | Could reject a solution error? |
|---:|---|---|---|
| 1 | Two short-coefficient forms; residual, `:99–110` | Only algebra of the face term; cannot test omitted tips | Yes, an algebra discrepancy; shared missing physics passes |
| 2 | Two steady-volume forms; residual, `:99–110` | No; endpoint volume identity | Yes, an algebra discrepancy; the displayed shared factor-of-two error passed |
| 3 | Hydraulic-J values required by three candidate tip factors, `:114–124` | No actual source integral or true tip factor was computed | No; a list of hypothetical values has no independent acceptance condition |
| 4 | Matrix determinant residual at three frequencies, `:126–137` | Operator algebra, not the tip matching factor | Yes, a matrix-algebra discrepancy; rerun arithmetic reproduces the residuals |
| 5 | Volume formulas after full-opening integration, `:126–141` | Endpoint normalization, not the short-time tip contribution | Yes; this comparison identifies the earlier volume factor of two |
| 6 | Nonzero matrix commutator, `:126–141` | No; diagnoses lack of constant simultaneous diagonalization | No as a solution check; no value of K or V is tested |
| 7 | Increasing grid resolution, `:145–157` | Yes, sampled initial tip intensity | Yes; numerical estimates could fail to stabilize, though the exact limit is not established |
| 8 | Decreasing diffusion length and extrapolation, `:145–157` | Yes, initial tip limit | Yes; estimates could fail to approach a consistent limit; the assumed fitting law remains unverified |
| 9 | Uncalibrated tip samples versus same-grid elastic calibration, `:145`, `:159–167` | Yes, tip-opening normalization | Yes; exposes the endpoint sampling bias |
| 10 | Linear-system relative residuals, `:148–155`, `:162–164` | Discrete equations only | Yes, a bad linear solve; cannot reject an incorrect operator or continuum extrapolation |
| 11 | Volume divided differences versus slope 1.2, `:157` | Yes, short-time volume growth | Yes; the supplied data contradict the claimed validation, as detailed above |
| 12 | Direct numerical reproduction output, `:169–187` | Prints the estimates rather than comparing them to an independent solution | No; wrong formulas or stored data can print consistently |
| 13 | Required keys present, `:191–198` | No | Delivery only; missing fields can fail |
| 14 | Printed JSON equals saved JSON, `:191–198` | No independent physical reference | Delivery only; divergent artifacts can fail |
| 15 | Numeric report entries, `:191–198` | No | Delivery only; invalid types can fail |
| 16 | General formula coverage, `:191–198` | Presence only, not mathematical correctness | Delivery only; absent strings can fail |

The fine-grid residuals and calibration are useful but do not establish an exact continuum intensity. No reported check detects the missing short-time tip contribution and correctly rejects the final coefficient; check 11 could have raised that issue but its evidence was misinterpreted.

## 6. Claims and reproduction

**Partly supported overall; all nine numerical reproduction claims are supported.** `AUTO_CHECKS.json` reports a successful rerun and all nine numbers printed. Its corrected checks reject only `V1_short_coeff`, identify the initial intensity as the tip-limit value, and accept both independent endpoint principles.

I copied the scripts and JSON to `scratch/reproduction/` and ran `python3 solve.py`: parsed stdout, regenerated JSON, and delivered JSON agree exactly, including the formula strings. The three runs of the delivered `tip_check.py` also reproduce their captured intensities and volumes to roundoff (maximum intensity difference `3.4e-16`, maximum volume difference `4.5e-16`). The slightly changed `~1e-12` linear residuals are ordinary rounding differences. The determinant residuals and commutator value were also independently reconstructed.

The unsupported part is the claimed slope validation in error 3. Also, `solve.py:10–15` does not freshly compute the continuum tip factor: it fits four stored extrapolated samples. The optional script reproduces only the three finest underlying runs, not all coarse runs or the estimated `ell=0.025` mesh correction. The answer accurately discloses this at `submission/ANSWER.md:30` and `:157`; it is not a mismatch or a claim of complete fresh reproduction. Arithmetic reproducibility does not validate the incorrect volume-growth physics. Numerical audit details are in `scratch/evidence.json`.

## 7. Questions and missing information

**Four listed entries** in the “Questions” section (`submission/ANSWER.md:23–26`), all interpretation statements rather than unanswered requests. The boundary/outer pressure distinction, full-jump volume normalization, and signed linear closure convention are sensible. Specifying contact or aperture is unnecessary for the given linear model, and the solver does not insist on it. The two-dimensional long-time/infinity caution is sensible and discussed below. None prevents the requested calculation.

## 8. Source naming

**None.** No submission file names a source paper or its authors. “Fourier-Laplace” and “Darcy flux” identify mathematical/constitutive concepts, not provenance for this crack problem. The authors and paper named in `KEY.md` are not credited to the solver.

## 9. Challenges and key assessment

The drained-tip versus undrained-outer distinction and warning about source-free elastic energy matching (`submission/ANSWER.md:80–93`) are consistent with the **corrected** key. The mixed operator and reproducible collocation support the distinction. They do not constitute a competing answer to the current key; the initial intensity is right against both the corrected key and the solver's tip interpretation. No independent defect in the corrected key was established in this audit.

One model/limit concern is explicitly raised: “A finite slit maintained at nonzero pressure in infinite two-dimensional space has no bounded strictly steady harmonic pressure field vanishing at infinity” (`submission/ANSWER.md:26`). The proposed resolution is to interpret final quantities as local limits of the transient solution, with spatial infinity and long time not commuting. The submission states this concern without supplying a proof or a global pressure-limit calculation. It is **checkable**: map the slit exterior to the exterior of a unit circle and solve the harmonic annulus problem `p_R(rho)=log(R/rho)/log(R)`. Taking `R` to infinity at fixed `rho` gives 1, rather than a globally decaying steady solution.

This is a sensible clarification, not an error or a developed alternative set of answers. Against the key, (b) is right; under the solver's local-transient interpretation, (b) is also right. A strictly steady global interpretation would require reconsidering existence, but the solver neither adopts it nor gives replacement K/V values. Accordingly there is no well-argued competing answer to mark contested; the checkable concern is retained in the profile below.

```json
{"answers": {"a": "partly", "b": "right", "c": "partly"},
 "supporting_state": "right",
 "method": "mixed Fourier-Laplace operator, periodic collocation and tip extrapolation", "method_class": "valid_alternative",
 "errors": {"decisive": 1, "side": 2, "presentation": 0},
 "checks_reported": 16, "checks_could_fail": 9,
 "claims_supported": "partly",
 "questions_listed": 4,
 "source_named": "none",
 "challenges": [{"claim": "No globally decaying strictly steady loading-2 pressure field exists in two dimensions; use the local transient limit", "checkable": true}]}
```
