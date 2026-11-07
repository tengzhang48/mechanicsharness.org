# Fixes after the audit of 25 October 2026

Audit: `audit_2026-10-25/AUDIT.md` (verdict READY AFTER FIXES). The statement in `task/` is kept exactly as the agents
received it; only the key and the checker were changed, and the saved runs were checked again. The files as run are in
commit `777e5ed`.

## Statement: kept, flaws named

- "K at the crack tips just after loading (t = 0⁺)" does not say which limit is taken first. The tip limit (r → 0 at
  fixed t > 0, then t → 0⁺) gives 0.807593, the key; the reverse order gives 1. The key grades 1 as partly right.
- The face shear stress is not prescribed in either loading; every reference assumed σ_xy = 0.
- `p_faces_0`, "immediately next to the crack faces at t = 0⁺", conflicts with the boundary value p = 0 on the faces; the
  intended value is the paper's Eq. 17.

## Key (`key/KEY.md`)

- The correction block (P374) set K1_0 = 0.807593 and V1_short_coeff = 2.440220, but the reference answer below it still
  gave the paper's values as the result. The reference answer now states the full solution: the drained tip matched to the
  undrained crack-scale field (K1_0), and the face term plus the tip term (A = 1.354055 + 1.086165).
- The paper's 0.75 (Eq. 29) and face-only 1.354 (Eq. 18) are described as its simplified solution and its face term, not
  as errors. The error was the first key, which took them as exact. An answer of 0.75 is recorded as the simplified
  solution.
- Eq. 27's prefactor written G√π/[√8(1 − ν)] (was √(πG)/…, dimensionally wrong).
- The exponent 0.743 is our reading for this material; Yu et al. (2020), Eq. 12, is for ν_u = ½, with α ≈ 0.735.

## Checker (`key/checks.py`)

- K1_0 band 1 % → 2 % and V1_short_coeff band 1 % → 5 %, with the label windows to match: a resolved FE route gives
  0.798–0.801 and 2.2–2.4 (FE_CHECK.md). The bands still exclude 0.75, 0.866, 1, 1.354 and 2/√π.
- Label `paper_heuristic_0.75` renamed `paper_simplified_0.75` (no saved run carried it).

## Saved runs checked again

`regrade_2026-10-26/SUMMARY.md`: no flag changes on the fifteen original runs. The three trials of 24 October (Kimi K3.1
high and max, Opus 5.5) had no saved automatic check; the fixed checker gives 5, 5 and 9 of 9 values, as the case page
reports. The blind grades and Grok reviews in `grading/` are as written at the time.

## Not applied

- The statement edits the audit proposes (the order of limits, σ_xy = 0, the definition of `p_faces_0`): the statement stays
  as run.
- Open before a public release: an independent check of the tip term 1.086165. Only the AI factorizations share it; our FE
  gives A within about 10 %.
- The audit text (25 Oct) calls the paper's values wrong under the corrected key; it was written before that wording was
  corrected and is kept as written.

## Corrections after GPT's review of the push (26 Oct)

- The labels `tip_limit_key` and `faces_plus_tips_key` now use the same relative bands as the pass checks (they used
  rounded absolute windows, 0.0162 and 0.122, slightly different from 2 % of 0.807593 and 5 % of 2.440220). No saved run
  changes.

## Decisions, 6 November 2026 (Teng)

- The open item above, an independent check of the tip term, is closed at the precision the check supports. A direct
  finite-crack solve that uses no factorization values gives 1.09 ± 0.01; a numerical Cauchy solve within the same tip
  model gives 1.0861646 ± 0.00002 (`independent_check_2026-10-29/REPORT.md`). The key keeps 1.086165. Its last digits rest
  on the shared tip-matching model, not on a second independent route. No grade depends on them: the term enters only
  through A = 2.440220, which is checked at 5 %.

## Checker, 6 November 2026 (found by Codex's review of the release)

- `key/checks.py` no longer takes JSON `true`/`false` as numbers: Python counts them as 1 and 0, so `K1_inf: true` and
  `K2_0: false` passed. A results file that is not a JSON object is now reported as unreadable (it used to crash the
  checker). Rerun on all eighteen saved answers: every output is identical to `regrade_2026-10-26/`. No grade changes.
