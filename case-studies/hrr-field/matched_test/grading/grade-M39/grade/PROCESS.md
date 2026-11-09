# Grading process

- 04:00:21 UTC — Read the task, problem, rubric, reference results, and submitted answer; began grading against the seven leading-field targets. Times are copied from `date -u +%H:%M:%S`.
- 04:00:31 UTC — Inspected all submitted Python files and the reference implementation, including the tensor reduction, flux calculation, shear convention, and three-principal-strain screening.
- 04:00:52 UTC — Copied code into scratch/ and reran solve.py, verify2.py, and solve_reference.py successfully. Confirmed the analytic outputs, independent J integration, displacement finite differences, and reference assertions.
- 04:01:17 UTC — Compared every reported output with the reference and saved relative errors in scratch/numerical_comparison.json. All seven pass. Recorded the script's alternative final interval and separated unsupported physical caveats from rubric requirements.
- 04:02:06 UTC — Wrote GRADE.md with the itemized rubric, 100/100 total, tolerance results, and scientific, support, and presentation notes. No reference correction was warranted.
- 04:02:22 UTC — Validated all 20 rubric rows sum to 100/100 and confirmed both required files and all seven passing tolerance comparisons. Grading complete; no web access or external process termination used.
