# Grading process

- 18:25:12 UTC — Started grading; read the task, mechanics problem, rubric, reference, and submission answer. Times are pasted from `date -u +%H:%M:%S`.
- 18:25:35 UTC — Reviewed all submitted code and the reference equations. Copied Python files into `scratch/` and launched the answer, face-residual, crack-limit, refinement, alternate-root, plane-stress, and reference checks there.
- 18:26:09 UTC — Reviewed completed check outputs and awaited the remaining alternate-route output. All seven initial checks then completed successfully. Confirmed p and R, integrator agreement, face residuals, the alternate-route result, and plane-stress comparison. Noted overstatements in the crack-limit and refinement precision claims.
- 18:26:52 UTC — Began a focused audit to measure the refinement difference at full precision and verify the tensile branch, energy power, and both requested tolerances.
- 18:27:07 UTC — Audit confirmed a 3.028688411177427×10^−12 refinement change, positive interior f, finite energy, and both quantities within tolerance. Applied the explicit partial-credit tiers: 10/20 for asserted dominance and 3/5 for scope; all other components receive full credit. Prepared GRADE.md with a total of 88/100.
- 18:28:11 UTC — Finished GRADE.md and verified its table sums to 88/100, both tolerance comparisons pass, copied input sources are unchanged, and rerun reference results equal the recorded JSON. Completed within the 25-minute limit. All started computations finished; no process was stopped or killed, and no web search or fetch was used.
