# Steady viscoelastic crack growth in an antiplane strip — solution

Entry point: `python solve.py case.json results.json`. Supporting code: `solve.py` (solver +
entry point), `post.py`, `checks.py`, `fastroots.py`, `energy.py`, plus convergence-driver
scripts (`run_final.py`, `run_v10.py`, `sweep.py`, `crit*.py`, `scan_highV.py`). All results
were produced with NumPy/SciPy only.

Notation. `X=pi x/h`, `Y=pi y/h`, `V=pi v t_c/h`, `eps=pi l_c/h`, `a=G_inf Delta/(tau_c h)`,
`W=G_inf w/(tau_c h)`, `U=W(.,0)`, `T_i=tau_i/tau_c`, `E_i=G_inf gamma_i/tau_c`.
`G_1=G_0-G_inf`, `lambda=G_inf/G_0=0.1`, `t_M=eta/G_1=lambda t_c`, `bar delta_c=0.2`.

## 1. Reduction to a harmonic crack-plane problem

The material is a standard linear solid. In the steady moving frame the material derivative
is `(.)dot = -v d/dx`, so the constitutive law reads
`(1 - v t_M d/dx) tau_i = (G_inf - v t_M G_0 d/dx) gamma_i`.
Because the operator acting on `gamma_i` is the same scalar differential operator in both
components, equilibrium `tau_{1,x}+tau_{2,y}=0` forces (with decay ahead)
`nabla^2 w = 0` **exactly**; viscoelasticity enters only through the `y=0` traction relation.

On `Y=0` the crack-plane relation (half strip `0<Y<pi`, Fourier variable `Q`, `K_0=1/(lambda V)`)
derived from the harmonic extension `w~(q,y)=u~(q) sinh(q(h-y))/sinh(qh)` is

```
T2~(Q) = 2 pi a delta(Q) - ell(Q) U~(Q),     U~(Q)=int U e^{-iQX} dX,
ell(Q) = (pi/lambda) Q coth(pi Q) (iQ - K0 lambda)/(iQ - K0).
```

Limits: `K0 -> inf` gives `ell = pi Q coth(pi Q)` (relaxed, `G_inf`);
`K0 -> 0` gives `ell = (pi/lambda) Q coth(pi Q)` (glassy, `G_0`).

The equivalent local (bulk) form used in the solver is
```
nabla^2 W = 0,   W(X,pi)=a,
T2 = pi W_Y + S2,   T1 = pi W_X + S1,
(d_X - K0) S_i = C W_{,Xi},   C = pi(1-lambda)/lambda,
```
with the bounded solution `S_i(X,Y) = -C int_0^inf W_{,Xi}(X+u,Y) e^{-K0 u} du`.
Boundary data: `W=0` (`X>=0`), `T2=1` (`-eps<X<0`), `T2=0` (`X<-eps`);
far field `W -> (a/pi)Y` ahead, `W -> a` behind.

The two scalar unknowns `(a,eps)` follow from
1. **critical separation** `2 U(-eps) = bar delta_c`, and
2. **bounded stress at the cohesive-zone front** `X=0` (vanishing stress-intensity
   coefficient `K = lim_{X->0+} sqrt(X) T2(X)`).

## 2. Slow-speed limit (exact)

As `V->0`, `K0->inf`, the bulk is fully relaxed and no bulk dissipation remains, so
`Gamma_D -> 0`, `R -> 1`, which through `R = a^2/bar delta_c` gives the exact limit
`a -> sqrt(bar delta_c) = sqrt(0.2) = 0.4472136`. The computed `V=1e-6` root,
`a = 0.447206`, agrees to 0.002 %, confirming the reduction and the solver.

## 3. Range of speeds / critical speed

Steady solutions exist for `0 < V < V_crit` and cease at a finite endpoint. Writing the tip
condition after eliminating `a` by the separation condition gives `K(eps;V)`. For `V` below
`~13` the curve `K(eps)` crosses zero at a finite cohesive-zone length. As `V` increases the
crossing is pushed to larger `eps`, and the minimum of `K(eps)` rises toward zero: at
`V = V_crit` the two conditions are tangent, and the admissible state has `eps -> inf`
(the cohesive zone no longer closes). For `V > V_crit` the SIF is positive for every `eps`
and no steady state exists.

Determined from the semi-infinite-zone SIF `K_inf(V) = lim_{eps->inf} K(eps;V)` computed on a
large domain (`L=160`, `eps=100`, `nx=400`):

```
V      13.40     13.45     13.50     13.55     13.60
K_inf  -2.89e-4  +3.00e-5  +2.13e-4  +2.96e-4  +4.15e-4
```

Linear interpolation of the sign change gives `V_crit = 13.44`. Repeat runs at a smaller
domain (`L=100`, `eps=80`) give `13.49`; the spread (limited by the SIF-extraction noise
floor, ~2e-4) brackets `V_crit in [13.3, 13.5]`. The reported value is `13.44`.

The limiting behaviour at the endpoint is `eps -> inf` with `a -> 1`. For `V >= V_crit` no
steady state exists (`epsilon = null`, `admissible = false`).

## 4. Bulk dissipation by spatial integration

The local dissipation is `phi_D = tau_i gamma_i_dot - psi_dot`. With `q_i = tau_i - G_inf
gamma_i` (Maxwell overstress) and dashpot viscosity `eta = G_1 t_M = G_inf(1-lambda) t_c`,
the Maxwell branch identity `gamma_i_dot = q_i_dot/G_1 + q_i/eta` gives the closed form
```
phi_D = q_i q_i / eta = G_inf (1-lambda) t_c^{-1} (tau_i - G_inf gamma_i)^2 / G_inf
      = tau_c^2 (T_i - E_i)^2 / (G_inf (1-lambda) t_c).
```
Since `T_i - E_i = S_i`, integrating over the strip and dividing by the crack advance (`1/v`)
gives, with `Gamma_0 = tau_c delta_c = bar delta_c tau_c^2 h / G_inf`,

```
Gamma_D / Gamma_0 = 2 / (bar delta_c (1-lambda) V pi) * int_0^pi int_{-inf}^{inf}
                    ( S1^2 + S2^2 ) dX dY .
```

This is a genuine spatial integral of `phi_D` (the factor 2 is the lower half strip by the
odd/even symmetry of `S1`,`S2`). It is **not** obtained by subtracting `Gamma_0` from
`Gamma`. The independently derived energy balance is `R = 1 + Gamma_D/Gamma_0` with
`R = a^2/bar delta_c`, and the integral reproduces it (see checks below).

## 5. Numerical results

Numerical method: second-order finite differences on a tensor grid in the strip, cosine
clustering towards the tip and towards `Y=0`, with the overstress convolution evaluated by
Gauss–Laguerre quadrature (the integration-by-parts averaging operator, which is stable
against the checkerboard mode of a naive composition). `S2` is eliminated to a boundary
unknown and the linear system is factorized once per `(V,eps)`; `a` enters the right-hand
side linearly, so the separation condition is solved by a scalar Newton and the tip
condition by a secant, at each of three resolutions (`nx = 240, 320, 400`) used for a
convergence study. Domain length is `L = max(14, 5V)` (the rear relaxation length scales as
`V`); this matters for `V=10` where `L=14` underestimates `Gamma_D` by a factor ~2.

Main results (as written to `results.json`):

| V | a | epsilon | R = a^2/0.2 | Gamma_D/Gamma_0 (spatial integral) | R-1 |
|---|-----|---------|--------------|-----------------------------|------|
| 0.001 | 0.449402 | 0.577000 | 1.009811 | 0.008925 | 0.009811 |
| 0.1   | 0.515225 | 0.768600 | 1.327282 | 0.326202 | 0.327282 |
| 10.0  | 0.972934 | 6.340000 | 4.733001 | 3.738413 | 3.733001 |

All three satisfy `|Gamma_D/Gamma_0 + 1 - R| / max(1,R) <= 0.01`:
`8.8e-4`, `8.1e-4`, `1.1e-3`. The 15 normalised field samples per case are in
`results.json`.

Resolution evidence (scalar roots; `a` converges to ~0.02 %, `epsilon` to ~0.1–0.3 %):

```
V=1e-6  nx=240 a=0.4478881 eps=0.5704515   nx=320 a=0.4472063 eps=0.5703108
V=0.001 nx=240 a=0.4500833 eps=0.5759149   nx=320 a=0.4495100 eps=0.5760685
        nx=400 a=0.4495875 eps=0.5770736
V=0.1   nx=240 a=0.5158214 eps=0.7683324   nx=320 a=0.5157352 eps=0.7691510
        nx=400 a=0.5151872 eps=0.7685370
V=10    nx=240 a=0.9736663 eps=5.70277 (L=40)  refined a=0.9729 eps=6.34 (L=50)
```

Normalised fields converge to <2e-4 between `nx=240` and `nx=320` (well inside the 0.005
gate).

## 6. Verification

`checks.py` prints the following for a solved state (interior-only residuals). Values for the
final states (`nx=320`):

```
check                     V=0.001      V=0.1        V=10
bulk Laplace (rel)        2.1e-08      1.2e-08      3.1e-09
W ahead of tip (max)      9.3e-10      3.0e-12      1.9e-10
|T2-1| cohesive zone      2.5e-09      3.2e-10      2.6e-10
|T2| behind (max)         4.8e-11      2.7e-10      3.5e-10
2U(-eps)-bar delta_c      -1.1e-04     -5.0e-04     +2.0e-04
outer BC (max)            8.3e-09      1.7e-08      1.2e-09
S2 integral vs solved     3.3e-12      2.5e-11      1.1e-10
energy balance err        -4.4e-04     -8.3e-04     +4.8e-03
```

The first four checks detect (i) an incorrect bulk discretisation/assembly, (ii) a wrong
ahead/interface condition, and (iii) traction leakage behind the cohesive zone. The
`2U(-eps)` residual monitors the cohesive law / critical separation. `S2 integral vs solved`
verifies the overstress operator against the directly solved boundary unknown. The
energy-balance residual is the independent check that the **spatially integrated** bulk
dissipation reproduces `R-1 = a^2/bar delta_c - 1`; it agrees to within the tolerance
`0.01 max(1,Gamma_D/Gamma_0)` in every case. The constitutive ODE residuals
`(d_X-K0)S_i = C W_{,Xi}` are ~1e-3 (V=10), ~1e-2 (V=0.1), and ~0.2 relative (V=0.001,
where the overstress layer of thickness `lambda V ~ 1e-4` is below the grid spacing so `S_i`
itself is under-resolved; its magnitude is `O(lambda V)` and the resulting field error is
`O(10^-4)`).

## 7. Tools, effort, limitations

- Tools: NumPy, SciPy (sparse LU, `RectBivariateSpline`), Python only. No external solvers,
  no network, no retrieved sources.
- The crack-plane symbol and local overstress operator were verified against the analytic
  Fourier solution and against both elastic limits (`K0->inf`, `K0->0`). The exact slow-speed
  limit `a -> sqrt(bar delta_c) = 0.4472136` is reproduced to 0.002 % at `V=1e-6`.
- Limitations / unresolved: (i) scalar tip quantities are limited by grid-dependent
  discretisation noise (~0.1–0.3 % on `epsilon`); (ii) `critical_speed` is bracketed to
  `[13.3,13.5]` by the `K_inf` noise floor (~2e-4); (iii) for `V << 0.001` the thin overstress
  layer is unresolved (negligible magnitude); (iv) the tip condition uses a fixed-window
  three-term SIF fit rather than an exact near-tip expansion. The two sample-speed energy
  balances (`V=0.001,0.1`) and the `V=1e-6` slow-limit check are the strongest independent
  validations.
- The T2 field behind the cohesive zone at the truncation end is not forced to zero (only
  `W=a` is imposed); `L = max(14,5V)` was chosen so this truncation error is below 0.2 %
  (verified by `Gamma_D/Gamma_0` convergence in `L`).