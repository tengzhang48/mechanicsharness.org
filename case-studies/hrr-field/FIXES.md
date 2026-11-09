# Fixes after the audit of 25 October 2026

Audit: `audit_2026-10-25/AUDIT.md` (verdicts T1 READY, T2 READY AFTER FIXES, T3 RESTATE). Only T1 and T2 are in this
folder; T3 invented specimens and measurements and is to be restated from the paper (draft
`audit_2026-10-25/T3_TASK_PROBLEM.md.proposed`, for Teng's approval). The statements in `task/T1` and `task/T2` are kept
exactly as the agents received them; only the rubrics were changed and a numeric checker was added. The files as run are
in commit `777e5ed`.

## Statements: kept, flaws named

- T1 line 30, "Cartesian differentiation of displacement must include rotation of the polar basis", warns about the
  pitfall that decides H and C. It helps without giving the answer.
- T2 is a variant: the paper applies its notch method only to a 90° notch at n = 3. The 120°, n = 4 answer is our
  computation, checked by two routes, by the paper's 90°, n = 3 exponent (0.22502 against 0.225) and by Fig. 9.
- Neither statement names the source, but every answer without the knowledge unit names "HRR".
- Neither statement asks for machine-readable numbers; the answers are free text.

## Keys

Unchanged. The audit reproduced T1's H and C to 2 × 10⁻⁹ and T2's p and R to 10⁻¹²; a scan of the whole finite-energy range
finds only one symmetric branch for T2, so "the leading local asymptote" is unambiguous.

## Rubrics (`key/T1/RUBRIC.md`, `key/T2/RUBRIC.md`)

- Both: whether all the numbers are at full credit is reported beside the total. R31 scored 55 with H and C wrong, and R23
  76 with both T2 numbers wrong.
- T2: a row for the governing equations (equilibrium and compatibility, 10 points), taken from the constitutive row
  (30 → 20). R23 used a wrong radial exponent in equilibrium and still received full constitutive points.
- T2: the 5 points for "the leading tensile branch" only if p is within the 0.003 band. R23 received them for a wrong
  eigenvalue.

## Checker (`key/check_P015.py`, new)

The audit's checker: it applies the rubrics' numeric bands to `{"T1": {p, H, C}, "T2": {p, R}}`. Every mutant the audit
tried fails it (neighbouring n, plane stress, C as πC, a missing n/(n + 1), the crack or elastic-wedge exponent for T2,
other notch angles).

On 8 November 2026 the checker was changed so that a value that is not a number (text, true/false) is reported as
"not a number", and a block that is not an object counts as empty, instead of stopping the check; its self-test output
on the saved runs and mutants is unchanged (`regrade_2026-10-26/check_P015_selftest.json`). An input without a T1 or T2
block now exits with status 1 ("nothing was checked") instead of printing an empty result (Codex's review, 8 Nov).

## Saved runs checked

`regrade_2026-10-26/SUMMARY.md`: with the numbers of each graded answer, T1 with-glm, T1 with-dsaz, T2 with-dsaz and T2
without-glm pass; T1 without-dsaz (H, C) and T2 without-dsaz (p, R) fail; two GLM runs timed out without an answer. These
agree with the graders' numeric rows. The totals are as graded.

## Record

`runs/P015-T1-without-dsaz/workspace/ANSWER.md` was overwritten seven minutes after the run's timeout by the run's own
helper session. The graded copy is now in `grading/grade-R31/submission/ANSWER.md`; `runs/P015-T1-without-dsaz/NOTE.md`
says so.

## Not applied

- The statement edits the audit proposes (a JSON output block for T1 and T2; optionally deleting T1's line 30): the
  statements stay as run.

## Corrections after GPT's review of the push (26 Oct)

- `key/check_P015.py --selftest` looked for its fixtures in `key/`; it now reads them from `audit_2026-10-25/`, and writes a
  file only when given `--out PATH`, so a rerun cannot overwrite the record in `regrade_2026-10-26/`. Its output is
  unchanged.
