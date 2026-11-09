# Grading process

- 04:01:13 UTC — Read the task, problem, rubric, reference results, and submitted answer. All work is local to this folder; no web access is used.
- 04:01:18 UTC — Inspected the submission and reference code. Began checking the constitutive reduction, numerical tolerances, screening conditions, and verification claims.
- 04:01:31 UTC — Copied the submitted solver into `scratch/` and reran it alongside the reference. The reference passed its assertions; the original submitted script stopped at an unavailable optional checker import.
- 04:02:00 UTC — Reran the submitted core calculations after removing only the unavailable optional checker import and its final check block in a scratch copy. No constitutive, field, inversion, screen, or contour calculation was changed. The numerical outputs and direct J checks reproduced the answer.
- 04:02:21 UTC — Completed the local audit: checked every reported number against its relative tolerance, all three principal strains, preservation of bulk compliance, direct contour values, and finite-r explanatory claims. Saved results in `scratch/audit_results.json`; no reference discrepancy was found.
- 04:03:04 UTC — Wrote and reviewed `GRADE.md`: 100/100 under the stated rubric, with separate lists of scientific errors, unsupported claims, and presentation problems. All seven requested numbers pass their tolerances.
- 04:03:19 UTC — Final consistency check confirmed that the 18 rubric rows sum to 100/100 and that all required report sections are present. Finished within the 25-minute limit.
