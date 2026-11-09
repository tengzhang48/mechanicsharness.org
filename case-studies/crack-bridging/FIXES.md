# Fixes after the audit of 25 October 2026

Audit: `audit_2026-10-25/AUDIT.md` (verdicts T1 READY, T2 READY AFTER FIXES, T3 RESTATE). Only T1 and T2 are in this
folder; T3 invented specimens and observations, handed over the decisive data and named the decisive issue, and is to be
restated from the paper (draft `audit_2026-10-25/T3_TASK_PROBLEM.md.proposed`, for Teng's approval). The statements in
`task/T1` and `task/T2` are kept exactly as the agents received them. Keys and rubrics are unchanged; a numeric checker was
added. The files as run are in commit `777e5ed`.

## Statements: kept, flaws named

- T1: none. Its numbers are parameters of a stated problem, inside the range the paper discusses
  (a/(δcE′/σb) = 0.302, Γ0/(σbδc) = 0.015).
- T2 is a variant of the paper's equal-arm splitting beam; the paper gives no solution for unequal arms, so the key is
  ours. The 20 mm starter cut does not enter the answer.
- T2's third output is δ(L/2), and δ(L/2)/δc depends only on Γ0/(σbδc), not on the arm stiffnesses, so it cannot detect the decisive
  unequal-arm error: averaged equal arms, or the upper arm alone, give the right δ(L/2). M and L do detect it.
- Neither statement asks for machine-readable numbers. The header the runs received is in each run's `TASK.md`; the
  current `setup_probe_runs.sh` added a rule after these runs were built.

## Keys

Unchanged. The audit recomputed T1's Σ and L by two routes (the paper's Eqs. 45–46 and an independent opening integral)
and T2's M, L and δ(L/2) by its own shooting and by energy, all to 10⁻¹³ or better. It also found a minor erratum in the
paper (p. 361: the equation numbers (35) and (36) are swapped), which affects no key.

## Checker (`key/check_P038.py`, new)

The audit's checker: it applies the rubrics' numeric bands (T1 ±0.3 %, T2 ±0.5 % for full credit) to
`{"T1": {Sigma_MPa, L_mm_per_end}, "T2": {M_Nm_each_strip, L_mm, delta_mid_um}}`.

On 8 November 2026 the checker was changed so that a value that is not a number (text, true/false) is reported as
"not a number" and a block that is not an object counts as empty, instead of stopping the check, and an input without a
T1 or T2 block exits with status 1 ("nothing was checked"), as for the HRR checkers after Codex's review; its self-test
output is unchanged (`regrade_2026-10-26/check_P038_selftest.json`).

## Saved runs checked

`regrade_2026-10-26/check_P038_selftest.json`: with the numbers of each graded answer, T1 R08, R24 and R30 get 50 of 50
and R03 (T1 without-dsaz) 0 of 50; all four T2 runs get 40 of 40. These agree with every grader's numeric rows. The
mutants all fail at least one item.

## Not applied

- The statement edits the audit proposes (a JSON output block; δ at x = 5.00 mm instead of δ(L/2), key 31.7700837 µm):
  the statements stay as run.

## Corrections after GPT's review of the push (26 Oct)

- `key/check_P038.py --selftest` looked for its fixtures in `key/`; it now reads them from `audit_2026-10-25/`, and writes a
  file only when given `--out PATH`, so a rerun cannot overwrite the record in `regrade_2026-10-26/`. Its output is
  unchanged.
- The note on δ(L/2) states the normalization: δ(L/2)/δc depends only on Γ0/(σbδc).
