# Grading process

- 03:58:30 UTC — Started grading; read TASK.md, TASK_PROBLEM.md, RUBRIC.md, REFERENCE.md, and submission/ANSWER.md. All work will remain in this folder; no web access will be used.
- 03:58:39 UTC — Read submission/solve.py, solve_reference.py, and reference_output.json; identified the rubric targets and reviewed the derivations.
- 03:58:47 UTC — Copied both scripts into scratch/ and reran them without modification. Both exited successfully; submission output is in scratch/submission_run.txt and the reference output is in scratch/reference_run.json.
- 03:59:25 UTC — Confirmed that all requested rounded values pass tolerance and that the numerical contour check is reproducible; reviewed written-equation errors separately against the scoring criteria.
- 03:59:55 UTC — Completed reproducible checks of all seven tolerances, code/reference agreement, direct J values, preserved bulk modulus, all three principal strains at the lower endpoint, and the erroneous prose calibration/energy expressions. Evidence is saved in scratch/grader_checks.json.
- 04:00:40 UTC — Wrote GRADE.md with every rubric subitem, the 100/100 total, scientific errors, unsupported claims, presentation problems, and the tolerance verdict. Ancillary prose defects do not change the correctly implemented scored calculations.
- 04:00:58 UTC — Final audit passed: all 20 rubric subitems sum to 100/100; all seven requested values pass their tolerances; required report sections are present; both original scripts match their unmodified scratch copies. Completed within the 25-minute limit, without web access or stopping any process.
