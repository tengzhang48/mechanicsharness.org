# Grading process

- 17:08:01 UTC — Read TASK.md, the mechanics problem, rubric, reference results and solution, and submission answer/code. No web access used.
- 17:08:11 UTC — Copied submission Python files to scratch/ for execution; checking numerical tolerances, implemented conditions, and verification claims.
- 17:08:30 UTC — Main and complex-exponential solvers reproduced the requested values. Numerical shooting reproduced the event profile, but its separate broad-bracket length search returned 48.705653 mm and emitted an fsolve convergence warning. Started the supplied continuation scripts and reference rerun.
- 17:09:01 UTC — Reference rerun exactly reproduced the stored numerical targets. Reviewed continuation outputs: correct first event, but post-failure traction and branch selection are not physically valid.
- 17:09:29 UTC — Completed quantitative audit against unrounded reference values, confirmed nonzero intrinsic front energy and rear pure-couple conditions, checked profile monotonicity and starter-cut invariance, and confirmed first-event shooting length with an admissible bracket. Awarded rubric credit from the submitted valid calculations, without counting failed or unimplemented checks.
- 17:10:29 UTC — Finished GRADE.md and validated all 13 scoring rows sum to 100/100, all three numerical tolerances pass, and all requested report sections exist. Retained rerun evidence in scratch/. Work completed within the 25-minute limit, with no web access and no process termination.
