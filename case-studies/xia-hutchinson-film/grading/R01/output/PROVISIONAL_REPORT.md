# PROVISIONAL REPORT

- UTC: 2026-10-15 (checkpoint_report)
- Role: independent scientific review of submission/ vs PROBLEM.md; reference used only to check source/model claims.

## Derived (reviewer, 1D wake)

- ell = sqrt(Ebar h / k), Ebar = E/(1-nu^2)
- gamma_A = tanh(H/(2 ell)); gamma_B = 2 tanh(H/(2 ell)) - tanh(H/ell)
- kI = sqrt((1-nu^2) gamma); isolated kI = 0.9539392014169457
- Independent FD of v''=v, v'(0)=v'(L)=-1: max |gamma_fd - tanh| ~ 7e-8 (independent_check.py / independent_check_output.json)

## Checks executed

1. Reviewer analytic+FD vs closed tanh: PASS (fd_err ~ 1e-8–7e-8)
2. Submitted 16 samples vs exact: n_agree at 1e-6 = 0; at 1e-3 = 3/16 (only gamma_B@1.5, gamma_B@2.0, kI_B@2.0). Typical absdiff ~ 0.002–0.005 (FEM bias, not 1e-6)
3. Coordinator replay (REPLAY_RECORD.json): exit 0; 11/16 display-rounding vs submitted JSON; 5 Case-B display mismatches; unguarded fem_energy2 import prints negative G / nan KI
4. Submission self-claims: Case B limits missing; KII argued not proven; J-integral and compute.py Finf rejected by solver

## Incomplete

- Final GRADE.md / grade.json / ANSWER.md / results.json / REPRODUCE.md not yet written
- No re-run of solver FEM (not required; coordinator already replayed)
- Source-recognition note (Xia–Hutchinson) to record as not claimed in submission
- Rubric items 1–7 judgments drafted in mind, not filed

## Uncertainty

- Submission uses 2D FEM energy differencing; far-field 1D reduction is the deciding math. Numbers are close in trend but fail 1e-6 and mostly fail 1e-3.
- Wrong far-field claim in ANSWER.md: “only bounded x1-independent equilibrium is u=0” — 1D cracked wakes are bounded nonzero u2(x2). That is a scientific error even though sampled gamma is nearby.
- Case B half-domain factor-of-2 and period accounting not independently re-derived from their FEM energy dump.
- Isolated-limit gamma_A→1 is correct; H/ell→0 they report gamma_A→0 linearly (consistent) but miss gamma_B ~ (H/ell)^3/4 and sequential shielding.
