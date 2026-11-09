| rubric item | points awarded / max | reason |
|---|---:|---|
| Finite geometry and coupling | 12 / 12 | Uses the coupled finite crack with c = a + L and closing traction on both a < abs(x) < c zones. The finite-crack tip and root-opening formulas match the reference; no small-zone approximation is used for the answers. |
| Intrinsic front work | 7 / 7 | Imposes K_tip = sqrt(E' Gamma_0) at each outer tip, separately from the connector closing traction. Neither doubles nor discards Gamma_0. |
| Elastic and local constitutive inputs | 6 / 6 | Correct plane-strain modulus E' = E/(1 − nu²) (3/3). The main formulation and supplied calculations use the constant 160 MPa closing traction and full relative opening with cutoff 9 micrometres (3/3). |
| Remote stress | 25 / 25 | Reported 133.73 MPa versus target 133.730286053811 MPa: absolute relative error 0.000213904%, within ±0.3%. The alternate rounding 133.730 MPa has the same error. |
| Zone length per end | 25 / 25 | Reported 4.086 mm per end versus target 4.08560913824553 mm: absolute relative error 0.009566793%, within ±0.3%. The alternate 4.0856 mm has error 0.000223669%, also within tolerance. |
| Verification | 15 / 15 | The supplied route B independently integrates the logarithmic opening kernel and agrees with analytical route A. Rerunning their unchanged A/B functions gives 133.73028605381106 MPa for both routes and L = 4.085609138245525 / 4.085609138245524 mm. Openings also agree at the four supplied branch locations to relative error below 1.3e−14. This satisfies the independent-elasticity check criterion; no credit is inferred from the unavailable spectral calculation or the incorrect supplementary limit check. |
| Conventions and event | 10 / 10 | Correct full-gap, per-end length, and single crack-plane-area conventions (5/5). Sets the opening at the original slit root equal to delta_c and identifies the branch followed from initiation; the supplied load sampling supports the monotonic branch before this event (5/5). |

**Total: 100 / 100.** The rubric's full-credit requirements are met despite the supplementary errors listed below; it specifies no additional deductions for those errors.

Execution evidence: the unmodified copied script fails at its unconditional `crack_bridging` import. To assess the mechanics actually supplied, only its existing constants, A/B function definitions, and A/B root statements were extracted unchanged into [scratch/extracted_routes.py](scratch/extracted_routes.py). No missing spectral implementation was supplied. Results are recorded in [scratch/audit_results.json](scratch/audit_results.json). Rerunning `reference.py` passed its assertions and reproduced the recorded key; no mechanical error in the reference was found. A grader calculation also confirms a positive load slope throughout the pre-event branch, with endpoint dimensionless slope 0.2150159321.

Scientific errors:

- Check 1's prose gives net kernel pressure as Sigma on the slit and −t on the connector zones. The zone pressure must be Sigma − t. Route B's code correctly superposes remote loading over the whole crack and closing traction on the zones, so the main calculation is unaffected.
- The code's short-zone opening check uses nonzero K0 and prints `delta_t_A(L)*Ep/(8*t*L)` with an asserted limit of 1. At its stated L = 0.01 mm it actually gives 1.705205746. For K0 = 0, the correctly normalized expression is `pi*delta*Ep/(8*t*L)` tending to 1 as L/a tends to zero. The finite-zone zero-K0 opening is logarithmic, so 8tL/(pi E') is an asymptotic expression, not an exact finite-L identity.
- For nonzero K0, 160 MPa is the large-c stress asymptote, not a strict ceiling. The submitted stress formula gives 160.69809 MPa at c = 100 mm. This does not invalidate monotonicity up to the requested failure event.
- The supplementary stress-intensity values 19.07 and 17.31 MPa sqrt(m) are inaccurately rounded: at the solution they are 19.08886453 and 17.32724635 MPa sqrt(m), respectively. Their difference still equals the prescribed intrinsic tip intensity.

Unsupported claims:

- The claimed spectral convergence sequence for N = 200 to 1600 and stress range 133.66 to 133.89 MPa is not demonstrated by the supplied files. The spectral implementation is absent, and the script contains no such refinement loop. Moreover, its printed route C stress is evaluated using the analytical `Sigma_of_c` formula, so that output alone would not establish an independent spectral stress result.

Presentation problems:

- The unconditional missing-module import prevents the submitted script from running its otherwise reproducible A/B calculations. The answer presents the checks collectively without distinguishing the supported A/B comparison from the unsupported spectral results and erroneous limit diagnostic.

Requested-number tolerance: **both Sigma = 133.73 MPa and L = 4.086 mm per end are within ±0.3% of the unrounded reference targets; every alternate answer rounding is also within tolerance.**
