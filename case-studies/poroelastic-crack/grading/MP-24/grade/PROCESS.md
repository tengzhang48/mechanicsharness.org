# Grading process

Times are pasted from `date -u +%H:%M:%S` (UTC).

- 19:29:34 — Read TASK.md, TASK_PROBLEM.md, KEY.md, AUTO_CHECKS.json, and all submission files. Confirmed corrected-key precedence, grading bands, and the requirement to assess challenges separately. No external sources accessed.
- 19:29:52 — Started `python3 solve.py` in `scratch/submission/`, with PID recorded in `scratch/submission/rerun.pid`. The original submission was not modified.
- 19:30:57 — Rerun output reproduced all nine numerical claims and both numerical checks. Inspected source at the cited lines to distinguish what the checks compute from what the answer infers.
- 19:33:14 — Completed `scratch/audit.py`: compared all values with the corrected key; checked the missing factor in the displayed volume calculation, the square-domain tail, and the offset slit grid. Derived and evaluated a finite-radius pressure-patch limit to test the loading-2 challenge: a locally uniform pressure need not imply zero total stress.
- 19:36:14 — Wrote GRADE.md with all nine requested sections, question-by-question judgments, cited error inventory, check scope and reproduction findings, and a separate assessment of the checkable loading-2 challenge. Applied the corrected key and its explicit partial-credit rule for the outer loading-1 intensity.
- 19:36:40 — Validated all nine sections, the required final JSON block and field values, error/check counts, and exact reproduction of the complete results object. Grading completed within the 25-minute limit; no process was killed and original submission files were left unchanged.
