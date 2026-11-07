# PROVISIONAL_REPORT — reviewer R04

UTC: 2026-10-15 (checkpoint_report)

## Derived / assessed so far

Submission is incomplete as a final package: no ANSWER.md, no solver results.json, no REPRODUCE.md. Mathematics is present in PROVISIONAL_REPORT.md, analytic_checks.py, numeric_eval.py, fd_check.py, run_all.py, numeric_rows.json.

Formulas in the partial write-up:
- ell = sqrt(h Ebar/k)
- gamma_A = tanh(H/(2 ell)), gamma_B = 2 tanh(H/(2 ell)) - tanh(H/ell)
- KI/(sigma0 sqrt(ell)) = sqrt((1-nu^2) gamma), KII=0
- sigma0c = sqrt(Ebar Gc/(ell gamma)); small-H and isolated-crack limits stated

Coordinator replay: python3 run_all.py, exit 0, 16 quantities within 1e-14 of original numeric_rows.json. Replay results.json is coordinator-generated, not original delivery.

## Checks executed (reviewer)

review_independent_check.py: independent trapezoidal quadrature of PROBLEM.md DeltaPi density on the 1D wake (not the solver tanh target). Agreement with tanh formulas ~1e-8 to 1e-10 (quadrature). Series: gamma_B ~ (H/ell)^3/4. 16/16 samples vs reference.json at 1e-12 (two kI_B diffs at 1e-16). Output: review_independent_check.py stdout.

Did not rerun fd_check.py in this review session; used submitted fd_output.txt and replay record. Solver STATUS.md had marked FD incomplete; later fd_check.py and fd_output.txt exist; replay ran run_all.py including FD.

## Incomplete

GRADE.md / grade.json / ANSWER.md / results.json / REPRODUCE.md. Full paper-source recognition scan of submission: no Xia/Hutchinson strings found. Item 6: distinguish solver FD vs reviewer quadrature vs coordinator replay. Item 7: Winkler-scope language in provisional report is thin.

## Uncertainty

Missing ANSWER.md is missing delivery, not repaired. FD path_spread ahead for Case A H=2 is 8.07e+00 in fd_output (J~0 ahead so relative spread noisy); |rel err| still <0.03. Command success ≠ physical film validation.
