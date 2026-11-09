# Grading process

- 17:04:37 UTC — Started grading; read TASK.md, TASK_PROBLEM.md, RUBRIC.md, REFERENCE.md, and submission/ANSWER.md.
- 17:04:45 UTC — Read submission/solve.py, reference.py, and the unrounded reference_results.json; checked the model and rubric subdivisions.
- 17:05:00 UTC — Copied submission/solve.py to scratch/solve.py and executed it successfully; saved output in scratch/solve_output.txt. Energy and statics checks pass; the finite-difference ODE residual is large.
- 17:05:14 UTC — Confirmed that the executed energy check earns full verification credit despite incorrect ancillary limiting-case statements.
- 17:05:56 UTC — Completed scratch/audit.py and saved scratch/audit_results.json: compared all reported and computed numbers with the unrounded key, checked the zero-intrinsic-work root, and diagnosed cancellation in the ODE stencil. All numerical rubric rows earn full credit.
- 17:06:42 UTC — Wrote GRADE.md and verified all 10 rubric rows, the 100/100 total, the three requested issue lists, and the tolerance statement. Grading completed within the 25-minute limit; all work remained inside this folder and no web access was used.
