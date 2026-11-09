# Grading process

- 17:04:40 UTC — Began grading; read the task, problem, rubric, reference, and submission answer. No web access used.
- 17:04:45 UTC — Read all four submitted scripts, the reference program, and unrounded reference results; created this log.
- 17:05:01 UTC — Reran unchanged copies of all submitted scripts and the reference inside scratch. Recorded outputs and exit statuses; confirmed reported main outputs and identified the main script's late exception and the separate alternative script's invalid opening profile.
- 17:05:48 UTC — Executed a grader-only audit against unrounded reference values. Verified that reversing the submission coordinate reproduces the reference profile and the positive intrinsic-front curvature; checked all three numerical tolerances and investigated additional scalar roots. Saved the evidence in scratch/audit_results.json.
- 17:06:29 UTC — Reviewed exact answer/code locations and applied the rubric's verification tiers. The submitted checks execute boundary residuals, but their energy check reuses the imposed resistance and the finite-difference ODE check prints a substantial residual. Grader-added checks do not earn submission verification credit.
- 17:07:28 UTC — Wrote GRADE.md with all rubric awards, reasons, the total, scientific errors, unsupported claims, presentation problems, and tolerance conclusion.
- 17:07:43 UTC — Validated all 11 rubric rows: awarded points sum to 94 and maxima sum to 100. Confirmed required lists, tolerance statement, process log, and saved execution evidence. Completed within the 25-minute limit.
