# Task: code the boundary checks in one anonymous answer (do not grade it)

Files: `TASK_PROBLEM.md` (what the solver was asked), `REFERENCE.md` (the reference solution and
reasoning, to judge whether a check is right), `submission/ANSWER.md`, and `submission/PROCESS.md`
(the solver's step log; lines that mention provided material are redacted).

For each of four **boundary checks**, decide whether the answer or the step log shows it:

- **B1 existence**: does the solver check that its solution form exists for these parameters
  (e.g. that a singularity or the assumed field exists, an eigenvalue is real and admissible,
  a finite-energy condition holds)?
- **B2 limiting cases**: does it check its result against a known limit through the same method
  (e.g. a linear-elastic or n = 1 limit, a closed-form special case)?
- **B3 applicability**: does it state the conditions under which its model holds and name physics
  the model leaves out (inner or outer limits, loading path, constraint, constitutive scope)?
- **B4 degenerate cases**: does it consider whether parameters sit near a change of structure
  (coinciding exponents, thresholds where a singularity appears or disappears)?

For each give: `status` = absent / present_wrong / present_right; `timing` = before / after /
unclear relative to the main numerical work (use the order of the step log; "after" includes
checks done only as final verification); `evidence` = at most two short quotes with file:line.
Judge only whether the check is present and physically right, not the rest of the answer.

Write `CODING.json`:
{"B1": {"status": ..., "timing": ..., "evidence": ...}, "B2": {...}, "B3": {...}, "B4": {...},
 "notes": "one or two sentences"}
and a short `CODING.md` with the same content in words. Keep `PROCESS.md` (times pasted from
`date -u +%H:%M:%S`). Rules: work only inside this folder; no web search or fetch. Time limit:
15 minutes.
Never stop or kill a process you did not start. Kill only by the PID you recorded when you
started it; never use pkill, killall, or pgrep-based loops (this machine is shared with other runs).
