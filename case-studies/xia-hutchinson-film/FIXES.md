# Fixes after the audit of 25 October 2026

Audit: `audit_2026-10-25/AUDIT.md` (verdict READY AFTER FIXES). The statement in `task/` is kept exactly as the agents
received it; a checker was added, and the saved runs were checked with it. The files as run are in commit `777e5ed`.

## Statement: kept, flaws named

- `PROBLEM.md` lines 3–4 carry coordinator text: "Draft for Teng's review. This is the proposed source-faithful
  reproduction task; source solutions and the proposed geometry extension are excluded." It tells the solver that the task
  reproduces a published source.
- The statement gives no accuracy for the sixteen values, while the record compared them at 10⁻⁶ absolute, and the method
  is open. S02's finite-element values (0.2–1.6 % off) counted 0 of 16.
- The source is recognizable, and its closed forms are printed in the paper (Eqs. 28, 30).

## Key

Unchanged. The audit recomputed every value by two routes (the closed forms, and a numerical solution of the strip
problem with quadrature of the stated functional, without tanh); they agree with the frozen reference to 2 × 10⁻¹⁶. The key
was written by GLM-5.3 and checked by Grok-4.6; the audit and the paper's printed equations now back it independently.

## Checker (`key/checks.py`, new)

- Reads the rows of `results.json` (or `numeric_rows.json`) and compares the sixteen values with
  `key/reference_frozen/reference.json` within 0.5 % relative; a missing value fails. Every physically wrong answer the audit
  tried fails at 0.5 %: no 1 − ν², no foundation energy, Case B per period, swapped spacings, ℓ/H read as H/ℓ, decay
  length 2ℓ. The 10⁻⁶ count is reported beside it.
- `key/RUBRIC.md` item 5 now names this count.

## Saved runs checked

`regrade_2026-10-26/SUMMARY.md`: S01, S04, D01 and K01 give 16 of 16 at both bands. S02 gives 7 of 16 at 0.5 % (0 at
10⁻⁶). S03 delivered nothing. The S01, S02 and S04 counts at 10⁻⁶ equal the record's `records/NUMERICAL_EVIDENCE.json`.
The Grok grades are as written at the time.

## Not applied

- The statement edits the audit proposes (delete the coordinator text, state an accuracy): the statement stays as run.

## Decisions, 6 November 2026 (Teng)

- **Release annotation** (drafted 29 October in `DECISION_PROBLEM_TEXT_2026-10-29.md`, approved 6 November). The solver
  input included the coordinator's draft/review wording at `task/PROBLEM.md` lines 3–4, quoted above. Those lines were
  part of the statement as run and remain unchanged. Treat these runs as source-reproduction trials with that cue; they
  are not evidence from a clean, uncued statement. This annotation changes neither the supplied problem nor the
  historical grades.
- The case page already says so where it shows the problem: "The task did tell every solver that it was a
  source-faithful reproduction and invited source recognition, so agents knew a published source existed. None named the
  paper; that does not rule out prior training exposure." The page is unchanged.
