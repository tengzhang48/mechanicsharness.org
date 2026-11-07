# PROVISIONAL_REPORT.md

UTC: 2026-10-16 (review checkpoint, minutes 8–10)
Phase: provisional reporting. No new investigations started.

## Completed judgments (provisional)

1. Wake formulation: **right** — ell=sqrt(h Ebar/k); A ahead u=0; A behind sinh strip; B ahead S=2H, behind S=H; total-traction faces.
2. Energy release: **right** on GA, GB, period/tip/h. Prose in ANSWER.md incorrectly says the linear prestress integral vanishes; compute.py still includes it and yields tanh.
3. Modes/KI: **right** — KII=0, kI=sqrt((1-nu^2) gamma).
4. Thresholds/limits: **right** — sigma_c, H->0 H/2ell and H^3/(4 ell^3), H->inf gamma->1.
5. Sixteen samples: **right** — 16/16 vs independent tanh/sqrt at 1e-12 vs reference.json.
6. Verification: **partly** — solver 1–8 share the 1D reduction; 2D J replay exit 0 with IntegrationWarning; this review has not yet executed a separate variational FEM.
7. Physical scope: **partly** — model limits stated; Winkler vs half-space / k-calibration not discussed.

## Independent checks this session

- Algebraic reduction from PROBLEM.md energy (not copied as a black box): gamma_A=tanh(H/(2ell)), gamma_B=2 tanh(H/(2ell))-tanh(H/ell).
- Direct comparison of submission/results.json rows to those formulas and reference/reference.json: 16/16.
- Coordinator factual replay: compute+checks exit 0; j2d/B_conv/A_H exit 0, IntegrationWarning present.
- Unexecuted: reviewer-owned P1 FEM of the full-width strip functional.

## Evidence

- submission/ANSWER.md, results.json, compute.py, checks.py, j2d.py
- factual_replay/CHECKPOINT_REPLAY_OUTPUT.txt, FINAL_REPLAY_check_j2d_cases_OUTPUT.txt
- reference/REFERENCE.md, reference.json

## Uncertainty / unfinished

- Independent FEM not run (budget checkpoint).
- Source recognition: solver says no external source; formulas match XH2000 §4.2 as frozen in reference/.
- Do not treat command success as physical validation.
