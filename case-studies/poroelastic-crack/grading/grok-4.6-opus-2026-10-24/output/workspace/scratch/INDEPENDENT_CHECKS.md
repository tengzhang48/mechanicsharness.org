# Independent checks (grader, not part of submission C)

Copied `C/code/crack_solver2.py` and saved numerical products into `scratch/` before any execution. No files under `C/` were written.

Scripts:
- `independent_checks.py` — closed forms, nine-value tolerances, re-extrapolation of saved large-s data, representative solver reruns
- `fit_and_identities.py` — polynomial-fit residuals, a(γ) bounds, 2-D capacity constant

Representative solver reruns (copied solver):
- s=100, N=20: bit-identical to `large_s_1e+02.json`
- s=1, N=24: Laplace-domain J identity relative residual 4.1e-10
- s=1e-4, N=16: matches `small_s.log` line 1

Did not rerun s=1e6–1e7 / N=150–190 or the Talbot campaign.

See script stdout from the grading session for the numerical tables.
