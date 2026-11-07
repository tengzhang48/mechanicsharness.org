# PROVISIONAL REPORT

- UTC: 2026-10-15 (checkpoint_report)

## Derived
Wake length ell=sqrt(h Ebar/k); strip fields from v''-v=0, v'(±W/2)=-1; U(W)=-tanh(W/(2ell)) in normalized units; GA, GB from translating-tip U_ahead-U_behind with N=1; Mode I KI=sqrt((1-nu^2) gamma); 16 samples at ell/H=0.5,1,1.5,2.

## Checks executed
- Coordinator replay: timeout 120s python check.py, exit 0 (factual, not a grade).
- Independent: python review_independent_check.py — solve_bvp + Simpson of stated energy; 16/16 agree at abs 1e-6 (diffs ~1e-12); closed_form_vs_bvp_max 1.098e-11.
- Solver CHECKS.txt: OVERALL PASS (eight self-checks); not treated as independent validation.

## Incomplete
None planned beyond reserved final packaging (GRADE/ANSWER already drafted).

## Uncertainty
No 2D tip-field K extraction. FD residual on BVP mesh is not a deciding check. Local-foundation vs real substrate not validated.
