# Grading process

- 03:40:37 UTC — Read the task, problem, rubric, reference solution, recorded reference output, and submission answer. Began checking the governing formulation and all seven signed numerical quantities against the 2% tolerance.
- 03:40:54 UTC — Finished reading all submitted scripts. Copied the scripts and reference code into scratch/ and started local reruns of final.py, verify2.py, and solve_reference.py.
- 03:41:08 UTC — All three reruns completed successfully. Checked every displayed answer value against the recorded reference and saved scratch/tolerance_check.csv. All seven pass; maximum relative error is 0.001044761%. Reproduced the flux, equilibrium, constitutive, and face-traction checks.
- 03:41:40 UTC — Confirmed the reference rerun exactly matches its recorded values and its two independent routes agree. Checked the first-pass solve2.py initialization and confirmed its erroneous negative-angle extrapolation; distinguished that draft and the earlier verify.py stress-scaling error from the corrected final implementations.
- 03:42:22 UTC — Wrote GRADE.md with itemized points, a 100/100 total, scientific errors in superseded drafts, unsupported claims, presentation issues, and the overall tolerance statement.
- 03:42:30 UTC — Validated the 17 rubric rows sum to 100 awarded and 100 possible points, confirmed seven passing numerical comparisons, and checked that all required report sections are present. Grading complete; all execution and generated files stayed inside this folder, with no web access.
