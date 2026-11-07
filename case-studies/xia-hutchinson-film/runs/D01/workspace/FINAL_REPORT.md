# FINAL REPORT — Periodic cracks in an elastic film on an elastic foundation

(Produced during continued analysis; ready for the minute-45-to-50 handoff.)

## 1. Model reduction and key results

Plane stress, thickness `h`, `Ebar=E/(1-nu^2)`, foundation `k`, reference
equibiaxial stress `sigma0`. Equilibrium gives

    (1-nu)/2 lap(u_a) + (1+nu)/2 d_a(div u) - (1/ell^2) u_a = 0,
    ell = sqrt(h Ebar / k) = sqrt(h E/(k(1-nu^2))).      (decay length)

Wake of a strip of width `W` (edges traction free):
`u1=0`, `u2 = B sinh((x2-W/2)/ell)`, `B = -sigma0 ell/(Ebar cosh(W/(2ell)))`,
`sigma22 = sigma0[1 - cosh((x2-W/2)/ell)/cosh(W/(2ell))]`, `sigma12=0`.
Strip-edge conditions: `u2' = -sigma0/Ebar`, `u1' = 0` (i.e. `sigma22=sigma12=0`).
Far ahead (both cases near the tips) and the Case-A far-ahead state is `u=0`;
Case-B far ahead is the old-crack wake (period `2H`).

Only **Mode I** is permitted (symmetry): `KII = 0`, `G = K_I^2/E`.

Energy accounting (per period, per unit advance, then /h per tip):

    G_A = (sigma0^2 ell/Ebar) tanh(H/(2 ell))
    G_B = (sigma0^2 ell/Ebar) [ 2 tanh(H/(2 ell)) - tanh(H/ell) ]

Normalized with `beta = H/(2ell)`, `nu=0.3`:

    gamma_A = tanh(beta),  gamma_B = 2 tanh(beta) - tanh(2 beta),
    kI_A = sqrt(0.91 tanh(beta)),  kI_B = sqrt(0.91 (2 tanh(beta)-tanh(2 beta)))

| ell/H | gamma_A | gamma_B | kI_A | kI_B |
|:-----:|:-------:|:-------:|:----:|:----:|
| 0.5 | 0.7615941560 | 0.5591607318 | 0.8324966558 | 0.7133276007 |
| 1.0 | 0.4621171573 | 0.1626401586 | 0.6484802334 | 0.3847109880 |
| 1.5 | 0.3215127375 | 0.0602425297 | 0.5409034952 | 0.2341382114 |
| 2.0 | 0.2449186624 | 0.0277201675 | 0.4720974293 | 0.1588249114 |

Critical reference stress for toughness `Gc`:
`sigma0_c(A)=sqrt(Gc Ebar/(ell tanh beta))`,
`sigma0_c(B)=sqrt(Gc Ebar/(ell(2 tanh beta - tanh 2beta)))`.
Limits: `H/ell->0`: `G_A ~ sigma0^2 H/(2Ebar)`, `G_B ~ sigma0^2 H^3/(4 ell^2 Ebar)`
(both `sigma0_c -> inf`); `H/ell->inf`: both `G -> sigma0^2 ell/Ebar`,
`sigma0_c -> sqrt(Gc Ebar/ell)`. Mechanical meaning: closely spaced cracks let
each tip relieve only a thin slice (Case B far less than A because the new
cracks only re-split an already relaxed wake); widely spaced cracks behave as
isolated semi-infinite cracks of the prestressed film-on-foundation.

Ordering: `G_A > G_B` for all finite `H`; `sigma0_c(B) > sigma0_c(A)`.

## 2. Verification actually executed (all PASS; 11 checks in results.json)

1D (same-reduction consistency): energy quadrature vs closed form (1.5e-16);
wake ODE residual (<7e-7); traction-free edges (0); independent numerical BVP
solve (<1e-10); integration-by-parts identity (1e-16); small-H limits
(`gamma_A~beta`, `gamma_B~2beta^3`); large-H limits (->1); `G=K^2/E` exact.

Independent 2D (different route: direct finite-difference solve of the PDE +
J-integral on a circular arc, no use of the wake energy formula):
- Case A `H=ell=1`: `2*J_arc = 0.46200..0.46226` for `r0=0.2..0.6` vs analytic
  `0.462117` (rel. err `~1e-4`), path independent.
- Case A H-dependence `H=0.5,1,2`: `G_num = 0.244833, 0.462170, 0.762159` vs
  analytic `0.244919, 0.462117, 0.761594` (rel. err `3.5e-4, 1.1e-4, 7.4e-4`),
  confirming the `tanh` shape.
- Case B: rel. err `1.9%, 0.94%, 0.62%` at `dx = 0.020, 0.010, 0.0067`
  (first-order in `dx`, ->0); Richardson extrapolation `G0 = 0.16271` vs
  analytic `0.162640`.

The 2D run shares only the mathematical model with the 1D work; it
independently confirms the coefficient, the per-tip/`h` normalization, and the
`H`-dependence.

## 3. Limitations / unresolved

- The model assumptions (plane stress, perfect bonding, fixed prestress during
  advance, linear elasticity) are the remaining limitation; no experimental
  `Gc` is supplied, so only the normalized scaling is reported.
- No numerical discrepancies were found.

## 4. Reproduction

    timeout 180 python3 compute.py && timeout 180 python3 checks.py
    timeout 300 python3 check_j2d_cases.py
    timeout 600 python3 check_B_conv.py
    timeout 500 python3 check_A_H.py

Machine-readable answer: `results.json`. Derivation: `ANSWER.md`.
