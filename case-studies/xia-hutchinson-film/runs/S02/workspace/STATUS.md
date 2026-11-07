# STATUS

- Time (UTC): 2026-10-15T13:48:00Z  (final reporting phase)
- Task: Solve PROBLEM.md (periodic cracks in elastic film on foundation).

## Completed (with saved evidence)
- Dispersion (derive3.py): branches s1^2=q^2+k/(hC), s2^2=q^2+2k/(hC(1-nu)),
  C=E/(1-nu^2). => ell = sqrt(Eh/(k(1-nu^2))) (branch-1 slowest).
- No x1-independent relaxed wake; disturbance is boundary layer width ~ell
  (derive8.py).
- CASE A FEM (fem_energy2.py, femA_H.py), grid-converged (dx=0.05 vs 0.075):
    ell/H=0.5: G_A=0.3481  gamma_A=0.7651  K_A/(s0 sqrt l)=0.8344
    ell/H=1.0: G_A=0.4235  gamma_A=0.4654  K_A/(s0 sqrt l)=0.6508
    ell/H=1.5: G_A=0.4423  gamma_A=0.3240  K_A/(s0 sqrt l)=0.5430
    ell/H=2.0: G_A=0.4494  gamma_A=0.2469  K_A/(s0 sqrt l)=0.4740
  (E=s0=h=1, H=1 units.)
- CASE A large-H limit: gamma_A -> 1.000, K_A/(s0 sqrt l) -> sqrt(1-nu^2)
  = 0.9539 as H/ell -> infinity (H=8,12,16 at ell=1).
- CASE A small-H trend: G_A -> 0; K_A/(s0 sqrt H) ~ 0.68 (rising) as H->0
  at ell=1; not fully converged.
- CASE B FEM (femb2.py, half-strip by symmetry), grid-converged (nx 600 vs
  1000 agree to ~0.3%):
    ell/H=0.5: G_B=0.2560  gamma_B=0.5642  K_B/(s0 sqrt l)=0.7158
    ell/H=1.0: G_B=0.1494  gamma_B=0.1642  K_B/(s0 sqrt l)=0.3865
    ell/H=1.5: G_B=0.0835  gamma_B=0.0612  K_B/(s0 sqrt l)=0.2359
    ell/H=2.0: G_B=0.0510  gamma_B=0.0280  K_B/(s0 sqrt l)=0.1598
  (factor 2 applied for half-domain.)
- Checks executed: FEM traction-free faces (sigma22~0 at wake face);
  crack-plane stress recovery to sigma0; energy vs crack length monotone;
  mesh convergence for A and B.

## Incomplete / unresolved
- H/ell->0 exact limit not fully converged (k=0 run interrupted).
- KII=0 only argued by symmetry, not proven.
- Case B large/small-H limits not computed.
- A discarded WRONG closed-form integral (compute.py) and an unreliable
  J-integral extractor were rejected; findings explicitly not used.
- results.json overwritten now with the FEM values above.

## Delivered
- ANSWER.md (final, with partial/labeled results).
- results.json (FEM-derived normalized values at sampled ell/H).
- PROVISIONAL_REPORT.md, PROVISIONAL + this STATUS, code, reproduction cmd.