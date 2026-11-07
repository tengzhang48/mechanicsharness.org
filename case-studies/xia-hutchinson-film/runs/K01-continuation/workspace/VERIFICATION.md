# VERIFICATION

Scope: provenance and meaning of every check, answering REVIEW_GAPS.md items 1–5.
No claim of expert validation or independent physical validation is made.

## 1. Check inventory with provenance

| Check | Kind | Prior attempt | This continuation session | Output |
|---|---|---|---|---|
| analytic_checks.py (SymPy, 8 checks) | internal algebra test | ran ~14:47 UTC, 8/8 PASS | rerun 19:41 UTC via run_all.py, 8/8 PASS | analytic_output.txt |
| numeric_eval.py (table + gamma_B identity) | internal consistency | ran ~14:48 UTC | rerun 19:41 UTC | numeric_output.txt, numeric_rows.json |
| fd_check.py (2D FD + Eshelby J) | independent numerical computation within the model | ran 15:06 UTC, PASS | independently rerun 19:41 UTC via run_all.py; values bit-identical to prior run | fd_output.txt |
| asymptotic-limit ratios (inline in run_all.py) | internal algebra/numeric | defined, replay reported by reviewer | executed 19:41 UTC, PASS | results.json checks block |
| fd_audit.py (path-spread diagnostic) | targeted diagnostic | did not exist | created and run 19:41 UTC (0.7 s) | fd_audit_output.txt |
| run_all.py (full replay) | reproduction command | reviewer: unchanged replay succeeded (artifacts not preserved) | executed 19:41 UTC, 8.73 s wall, ALL CHECKS PASS, results.json regenerated | replay_output.txt, results.json |

"Independent" is used narrowly: the FD route is independent of the closed-form algebra
(own grid, own solver, own J extraction), but NOT independent of the model equations or of
the analytic far states (see §3).

## 2. Review-gap item 1: deliverables and the stale FD status

- ANSWER.md, results.json, REPRODUCE.md now exist (this folder).
- The prior STATUS.md (15:00 UTC) listed the FD work as incomplete; fd_output.txt
  (15:06 UTC) and run_all.py (15:07 UTC) were produced AFTER that checkpoint, so the FD
  check did complete in the prior attempt and the STATUS line was simply stale. Confirmed
  in this session: the rerun reproduces the saved fd_output.txt values exactly
  (G_fd=0.420499/0.692932/0.222872/0.147913/0.508152; ratio 3.93).

## 3. Review-gap items 2 & 3: what the FD agreement tests and what it cannot validate

The FD calculation prescribes the analytical far states as Dirichlet data at the domain
ends and compares the extracted G with the closed-form tanh target. Therefore:

It DOES test:
- the translating-tip energy bookkeeping: that G per advancing tip equals
  (W_ahead - W_behind)/h with the stated period widths (H for A, 2H for B) and tip counts
  (one per period in both cases), including the h normalization;
- that the closed-form G coincides with the Eshelby J of an actual 2D field solution of
  the model BVP (crack-face and symmetry edge conditions enforced by ghost-node
  elimination), to <= 0.13% over 5 configurations, with grid-convergence ratio 3.93
  (consistent with the second-order scheme) and path-independent J (see §4);
- the J-flux identity J_behind = -G h/2 (match to 6e-5 absolute).

It CANNOT validate:
- the strip-ODE far states themselves: they enter as prescribed boundary data, and the
  comparison target is built from the same states — agreement there is by construction;
- the model PDE: the FD solves the same equations the analytics reduce, so an error in the
  model statement (local foundation law, plane-stress reduction, prestress bookkeeping
  inside DeltaPi) would pass undetected by this route;
- anything about a physical specimen: no material data, no calibration, no experiment.

Shared assumptions between verification routes: the model equations, the far-field
bounded equilibria, and G = KI^2/E for the KI normalization. The SymPy and limit checks
share, in addition, the closed forms themselves (they are internal algebra tests).

## 4. Review-gap item 4: path-spread 8.07 audit (fd_audit.py, fd_audit_output.txt)

Case A, H = 2.0, ahead region. The statistic is (max-min)/|mean J| over the ahead window.
Measured this session (single re-solve, 0.7 s):
- J_ahead mean = -9.81e-10 (continuum target: exactly 0 in the undisturbed state);
- absolute spread = 7.92e-9, giving 7.92e-9 / 9.81e-10 = 8.07 — a divide-by-near-zero
  artifact of the relative normalization, not a path-independence failure;
- as an absolute error relative to the released energy (G h/2 = 0.3465): 2.3e-8;
- the behind region shows relative spread 1.6e-6 (absolute 5.6e-7) — the ahead region is
  in fact better converged in absolute terms.
Conclusion: path independence holds to O(1e-6) absolute everywhere; the appropriate error
measure near J = 0 is the absolute spread, reported above.

## 5. Review-gap item 5: physical scope

- The local Winkler foundation (f = -k u) and the plane-stress membrane reduction are the
  GIVEN model. Substrate stiffness k has not been calibrated for any real specimen; no
  dimensional material data were supplied or used (all reported values are normalized).
- The dense-spacing asymptotics (sigma0c ~ H^{-1/2} Case A, ~ H^{-3/2} Case B as
  H/ell -> 0) are asymptotic properties of this reduced model at fixed k. They are not
  automatically a validation of, or prediction for, a real substrate: real compliant
  substrates have finite thickness, non-local stiffness, and possible debonding or
  inelasticity, all outside this model.
- Supported physical conclusions are therefore conditional on the model: the length
  ell = sqrt(h Ebar/k), the closed-form G_A(H), G_B(H), pure mode-I symmetry, and the
  threshold formula sigma0c = sqrt(Ebar Gc/(ell gamma)). Anything beyond — substrate
  calibration, crack-path stability, dynamics, large-deformation or debonding regimes —
  requires additional information or validation not present in this task.

## 6. Declared tolerances and measured outcomes (this session's rerun)

- FD |rel err| < 3%: measured max 0.13% (5 configs) — PASS.
- Grid convergence ratio consistent with 2nd order (~4): measured 3.93 — PASS.
- SymPy residuals: exactly 0 (8/8) — PASS.
- gamma_B two forms: |diff| < 1e-15 — PASS.
- Limits: gamma_A/x and gamma_B/(2x^3) -> 1 within 1e-4 at x = 1e-3 (measured 3.3e-7 and
  2e-6); gamma(x=20) = 1 to machine precision — PASS.
- Full replay wall time: 8.73 s (within the 120 s per-computation guidance).
