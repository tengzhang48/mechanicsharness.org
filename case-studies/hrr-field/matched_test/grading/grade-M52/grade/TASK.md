# Grading instructions

You grade one anonymous submission to a mechanics problem. The folder contains:
`TASK_PROBLEM.md` (what the solver was asked), `RUBRIC.md` (100 points), `REFERENCE.md` with the
reference code and its recorded results (the key), and `submission/` (the solver's `ANSWER.md` and
any code it wrote).

1. Grade strictly against `RUBRIC.md`, using the reference as the key. Award partial credit exactly
   as the rubric states. Check every requested number against its stated tolerance.
2. Where the submission's numbers depend on its code, you may rerun it (copy it into `./scratch/`).
   Credit what the answer and its code actually show, not what they claim.
3. If you believe the reference itself is wrong somewhere, say so with evidence, and grade both
   ways (against the reference, and against your corrected view).
4. Some submissions had extra material available and some did not. Ignore any mention of it; grade
   only the mechanics in the answer.

Write `GRADE.md`: a table `| rubric item | points awarded / max | reason |`, the total out of 100,
then lists of (a) scientific errors, (b) unsupported claims, (c) presentation problems, and one line
on whether the requested numbers are within tolerance. Keep `PROCESS.md` (times pasted from
`date -u +%H:%M:%S`). Rules: work only inside this folder; no web search or fetch. Time limit: 25
minutes.

Never stop or kill a process you did not start. Kill only by the PID you recorded when you started it; never use pkill, killall, or pgrep-based loops (this machine is shared with other runs).
