| rubric item | points awarded / max | reason |
|---|---:|---|
| Finite geometry and coupling | 12 / 12 | Uses the coupled finite crack with half-length c = a + L and closing traction on both connector zones. The finite-crack intensity and root-opening formulas match the reference; no small-zone approximation is imposed. |
| Intrinsic front work | 7 / 7 | Enforces K_tip = sqrt(E′Γ₀) = 1.76161818 MPa√m at each outer front, with Γ₀ = 22 J/m² separate from connector work and counted once. |
| Elastic and local constitutive inputs | 6 / 6 | 3/3 for E′ = E/(1−ν²) = 141.0590278 GPa; 3/3 for constant closing traction 160 MPa acting against the full gap, with cutoff 9 µm. The opening normalization is correct. |
| Remote stress | 25 / 25 | Reported 133.7 MPa differs from 133.730286053811 MPa by 0.0226471%; the parenthetical 133.73 MPa differs by 0.000213904%. Both satisfy ±0.3%. Rerun code gives 133.73028605381106 MPa. |
| Zone length per end | 25 / 25 | Reported 4.09 mm per end differs from 4.08560913824553 mm by 0.1074714%; the parenthetical 4.086 mm differs by 0.00956679%. Both satisfy ±0.3%. Rerun code gives 4.085609138245525 mm per end. |
| Verification | 15 / 15 | Executes and interprets an independent dislocation-density elasticity integration: the zero-tip-work opening agrees to 5.80×10⁻¹¹ relative error, and the total root opening is 8.999999999600 µm. The separate spectral calculation also executes and gives 133.69987475 MPa and 4.07875639 mm per end. These meaningful checks earn full credit despite the isolated erroneous initiation-limit sentence below. |
| Conventions and event | 10 / 10 | 5/5 for full-gap, per-end length, and crack-plane area conventions. 5/5 for identifying δ(a⁺) = δ_c as first failure: the increasing root opening gives a unique cutoff, the computed profile decreases toward the outer tip, and the submitted branch table shows increasing remote stress from initiation through this event. |

**Total: 100 / 100.**

Scientific errors:

- The final check in ANSWER.md incorrectly states that K_tip tends to zero as c tends to a. Along the branch actually implemented, K_tip = K₀ = 1.76161818 MPa√m remains nonzero; the zone-root opening tends to zero. The earlier tip condition and initiation stress are correct, so this sentence does not change the calculated event or negate the successful independent checks.
- The spectral stress agreement is stated as “within 0.02%”; the rerun discrepancy is 0.0227408%. Its length discrepancy is 0.167729%, consistent with the stated 0.17%. This small overstatement does not affect either requested answer's numerical credit.

Unsupported claims:

- The refinement results bracket the exact tip intensity and opening evaluated at the closed-form solution. Bracketing of the independently solved stress and length across resolutions is not demonstrated: final_check.py solves that system only at M = 4096.

Presentation problems:

- solve.py labels a profile sample “along zone,” although it samples the original slit from the center to just inside x = a. final_check.py supplies the full profile through the connector zone.
- The “pure-remote” anchor in solve.py prints a label without executing a comparison; it earns no verification credit. Units and the per-end convention in the answer are otherwise clear.

The reference rerun passes its assertions and reproduces both recorded targets exactly; no scientific error in the reference was found. Execution evidence is saved in [scratch/submission_run.txt](scratch/submission_run.txt), [scratch/audit_results.txt](scratch/audit_results.txt), and [scratch/reference_run.json](scratch/reference_run.json).

**Requested-number tolerance: both Σ and L are within ±0.3% of the unrounded reference targets, including every displayed rounding of the requested results.**
