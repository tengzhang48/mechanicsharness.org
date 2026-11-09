# Grading process

- 18:25:12 UTC — Started grading; read `TASK.md`, the problem, rubric, reference, and submission answer. No external sources will be used.
- 18:25:23 UTC — Read all submitted Python files and the recorded reference checks. Identified `-np` in the submission's equilibrium equations, inconsistent with its stress ansatz `r^-p`; computed absolute errors in the two requested numbers.
- 18:26:08 UTC — Copied code into `scratch/`. Direct checks reproduced the linear-notch exponent and the submission's second, infinite-energy mode. Background launch attempts produced no results, so reruns use retained command sessions and record their own PIDs.
- 18:26:32 UTC — The reference rerun exactly reproduced the supplied result JSON. The submission's n-continuation reproduced its reported exponent and ratio, but returned mesh-limit failures for n > 1.
- 18:27:44 UTC — Completed the submission's n = 1 angle-continuation rerun; its final near-crack field failed with residual about 109. Three grid-refinement rows reproduced the reported digits while failing their requested convergence tolerances.
- 18:28:21 UTC — All four original `conv.py` refinement runs completed; all returned status = 1, with achieved residuals 4.8e-9 to 6.2e-9 and identical reported digits. Applied the rubric's explicit partial-credit levels and prepared `GRADE.md` with total 76/100. No processes were killed.
- 18:28:54 UTC — Verified the six rubric rows sum to 76/100 and the required error lists, tolerance verdict, and rerun evidence are present. Completed grading within the 25-minute limit, using only this folder and no web access.
