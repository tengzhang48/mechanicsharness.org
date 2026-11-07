# PROVISIONAL REPORT (minute ~4-5, UTC 2026-10-15T21:22:33Z)

Status: analysis essentially complete for the analytic core; deliverables being
finalized. Number below are provisional but already cross-checked.

## Completed
- Decay length: ell = sqrt(h Ebar / k), Ebar = E/(1-nu^2).
- Only Mode I is permitted by the symmetry of both geometries; KII = 0.
- Far-behind wake (strip of width W, edges x2=const traction free):
  u1 = 0, u2 = B sinh((x2-W/2)/ell), B = -sigma0 ell/(Ebar cosh(W/(2ell))),
  sigma22 = sigma0[1 - cosh((x2-W/2)/ell)/cosh(W/(2ell))], sigma12 = 0.
  Strip-edge conditions: u2' = -sigma0/Ebar, u1' = 0, i.e. sigma22 = sigma12 = 0.
- Energy per period per unit advance (release), then per tip / h:
  GA = (sigma0^2 ell/Ebar) tanh(H/(2ell))
  GB = (sigma0^2 ell/Ebar) [2 tanh(H/(2ell)) - tanh(H/ell)]
- Normalized: gamma_A = tanh(beta), gamma_B = 2 tanh(beta) - tanh(2 beta),
  beta = H/(2ell); kI = sqrt(gamma (1-nu^2)).
- Limits: GA ~ sigma0^2 H/(2Ebar) (H->0), -> sigma0^2 ell/Ebar (H->inf);
  GB ~ sigma0^2 H^3/(4 ell^2 Ebar) (H->0), -> sigma0^2 ell/Ebar (H->inf).

## Provisional sampling table (nu = 0.3)
| ell/H | gamma_A | gamma_B | kI_A | kI_B |
|------|---------|---------|------|------|
| 0.5 | 0.761594 | 0.559161 | 0.832497 | 0.713328 |
| 1.0 | 0.462117 | 0.162640 | 0.648480 | 0.384711 |
| 1.5 | 0.321513 | 0.060243 | 0.540903 | 0.234138 |
| 2.0 | 0.244919 | 0.027720 | 0.472097 | 0.158825 |

## Checks executed (all PASS so far)
- direct energy quadrature vs closed form: max rel err 1.5e-16
- wake ODE residual: < 7e-7
- traction-free edges: 0
- independent numerical BVP solve vs closed form: < 1e-10
- integration-by-parts identity: rel 1e-16
- small/large H limits: match expansions
- G = K^2/E consistency: exact

## Unresolved / pending
- Genuinely independent 2D verification of the energy-release coefficient
  (J-integral on a 2D numerical solve) not yet executed; the checks above
  share the 1D reduction, so they do not independently confirm the factor.
- Final ANSWER.md / results.json / REPRODUCE.md have now been saved (this
  checkpoint). All requested quantities are populated; the only outstanding
  item is the independent 2D cross-check noted above (explicitly flagged in
  ANSWER.md §8 and results.json "unresolved").
