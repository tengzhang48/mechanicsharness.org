# PROVISIONAL REPORT — continuation checkpoint (minute 5 handoff)

- UTC time: 2026-10-15 19:39:46Z (continuation session started 19:35:01 UTC; this is the
  reserved minute-5-to-6 reporting phase).
- Phase: checkpoint_report. No new calculation was started in this window; everything
  below is saved work.

## Deliverables already saved this session
- STATUS.md and CONTINUATION_CHECKPOINT.md (saved first, 19:35 UTC).
- ANSWER.md — full derivation from the existing work: ell = sqrt(h Ebar/k) with the two
  eigen-branches; far-ahead/far-behind wake states for Cases A and B with all strip-edge
  conditions; translating-tip/J energy argument with period width, tip count and h
  normalization; G_A = (sigma0^2 ell/Ebar) tanh(H/2ell), G_B = (sigma0^2 ell/Ebar)
  [2 tanh(H/2ell) - tanh(H/ell)] = 2t^3/(1+t^2) form; pure mode-I symmetry, KII = 0,
  KI/(sigma0 sqrt(ell)) = sqrt((1-nu^2) gamma); sigma0c = sqrt(Ebar Gc/(ell gamma));
  H/ell -> 0 and -> inf limits with mechanical meaning; numerical table; check inventory
  with what each check can and cannot detect; model limitations; open questions.
- results.json — normalized quantities, sampling locations, critical-stress formulas,
  check list with honest status (prior-executed vs replay-pending), model limitations.
  Note: a later `python3 run_all.py` replay will overwrite results.json with fresh
  executed-now outputs; content expected identical for the analytic rows.
- REPRODUCE.md — single command `python3 run_all.py`, step list, tolerances, requirements,
  provenance and shared-assumption notes.

## Checks actually executed so far (provenance)
- Prior attempt (saved files, executed 2026-10-15 14:47–15:07 UTC):
  - analytic_checks.py: 8/8 PASS (analytic_output.txt).
  - numeric_eval.py: table + numeric_rows.json; two-form gamma_B assert < 1e-15.
  - fd_check.py: COMPLETED at 15:06 UTC (fd_output.txt) — 5 configurations |rel err|
    <= 0.13% (tolerance 3%), convergence ratio 3.93, J_behind flux match 6e-5. The old
    STATUS.md (15:00 UTC) calling FD "incomplete" was stale: the run finished after that
    checkpoint. Confirmed top-level files identical to previous_attempt/ copies.
  - run_all.py written 15:07 UTC; reviewer reports the unchanged replay succeeded; no
    results.json was preserved from the prior session.
- This session (as of 19:39:46 UTC): NO new computation executed yet. Assembly of
  ANSWER.md / results.json / REPRODUCE.md from saved outputs only.

## Remaining work / gaps
- [ ] Independent rerun in this session: `python3 run_all.py` (includes FD; the one
      expensive step, est. tens of seconds to ~2 min) — audit/reproduce the numerical check.
- [ ] Audit of path_spread ahead = 8.07 (Case A, H=2.0): relative spread normalized by
      |mean J| where J ~ 0 ahead; report absolute spread (~1.6e-6, i.e. ~5e-6 of the
      released energy) — cheap diagnostic, may reuse saved output or one small re-solve.
- [ ] VERIFICATION.md — prior vs new checks, what FD agreement tests and cannot validate,
      shared assumptions.
- [ ] Final STATUS.md and final reporting (minutes 15–18, 19:50–19:53 UTC).

## Risk notes
- If the FD rerun exceeds the 120 s computation guidance, the saved prior fd_output.txt
  remains the evidence of record and the rerun will be reported as not completed.
- All REVIEW_GAPS.md items are mapped: (1) final deliverables now exist; FD-status
  discrepancy explained above; (2) rerun pending; (3) FD scope explained in ANSWER.md §5
  and to be detailed in VERIFICATION.md; (4) path-spread audit pending; (5) model-scope
  statement in ANSWER.md §6.
