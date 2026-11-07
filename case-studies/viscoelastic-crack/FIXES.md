# Fixes after the audit of 25 October 2026

Audit: `audit_2026-10-25/AUDIT.md` (verdict READY AFTER FIXES). The statement in `task/` (`PROBLEM.md`, `CONTRACT.md`,
`case.json`) is kept exactly as the agents received it; only the scorer was changed, and the saved submissions were scored
again. The files as run are in commit `777e5ed`.

## Statement: kept

The audit found no ambiguity, erratum or invented narrative; every value in the key was recomputed independently. The
statement is recognizable as Hui and Long (2026), and `CONTRACT.md` asks solvers to record recognition.

## Scorer (`key/score_submission.py`)

`CONTRACT.md` announces the accuracy gates to the solver, so those gates are kept as announced and still decide
`all_required_numeric_gates_pass`. Two kinds of check are added and reported beside them, as `audit_checks_2026_10_26`:

- **Bulk dissipation in a relative band.** The announced band, 0.01 max(1, Γ_D/Γ₀), passes Γ_D = 0 at V = 0.001, where
  the reference is 0.0073; P01 and P04 passed there with errors of +22 % and +49 %. The added check uses 5 % of the
  reference or 2 × 10⁻⁴, whichever is larger, and the energy balance |Γ_D/Γ₀ + 1 − R| within 5 % of R − 1 plus 2 × 10⁻⁴.
- **The submitter's own numbers against the source's relations,** a = (2/π) acos(e^(−ε/2)) (Eq. 18) and R = a²/δ̄ (Eq. 22),
  each within 0.5 %. P01's inconsistency in a (2.45 %) had been found only by hand.

Checked: the reference passes; a submission with Γ_D = 0 at V = 0.001 passes the announced gates and fails the added ones.

## Saved submissions scored again

`regrade_2026-10-26/`: the announced gates give the same result as the saved reports for P01 (8 failed checks) and P04
(31 failed). P03 (GPT-6.1 Sol) passes the announced gates and the added checks. The added checks fail for P01 (bulk and
energy at V = 0.001, a from ε at V = 0.001 and 0.1) and for P04 (bulk and energy at V = 0.001 and 10; at V = 10, which it
declared inadmissible, it gave no ε, a or R). P02 and P05 (Kimi) delivered no results.

## Rubrics

Two rubrics exist: `key/SCORING.md` (100 points), frozen before the runs, and `key/SCORING_ANALYTICAL_V2.md` (70/20/10),
added after launch and after GPT's delivery at the user's request. SCORING.md is the rubric of record; V2 is a post-hoc
reading. Grok's review gives GPT 100/100 under both (`records/GROK_GRADING_STATUS.md`).

## Not applied

- The audit's change to the band announced in `CONTRACT.md`: the contract stays as run.
- Absolute paths of this machine remain in the harness scripts; the repository is private.

## Corrections after GPT's review of the push (26 Oct)

- The added relation checks evaluated exp(−ε/2) and acos(…) before testing ε, so a negative ε (−0.5, −2000) raised an
  error instead of a report. They now require positive ε, a and R first, as the contract asks, and report a failure
  otherwise. The saved submissions score as before.
