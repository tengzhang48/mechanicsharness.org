# ANSWER.md — Periodic cracks in an elastic film on an elastic foundation

Assumptions and units: plane stress, isotropic, Young's modulus `E`, Poisson
ratio `nu`, film thickness `h`, foundation stiffness `k` [Pa/m], reference
equibiaxial stress `sigma0` [Pa]. `Ebar = E/(1-nu^2)`. Numerical values use
`nu = 0.3`. Lengths in metres; `ell` (defined below) is the only length scale
of the wake. All quantities are reported in the problem's requested
normalization.

The equations, geometry and boundary conditions are taken exactly from
PROBLEM.md. No external source was consulted.

---

## 1. Governing equation and decay length

With `sigma = sigma0 I + Delta sigma` and
`Delta sigma_{ab} = Ebar[(1-nu) eps_{ab} + nu eps_{gg} delta_{ab}]`,
`eps_{ab} = (u_{a,b}+u_{b,a})/2`, equilibrium `h sigma_{ab,b} - k u_a = 0`
(the constant `sigma0 I` has zero divergence) becomes

    (1-nu)/2 * lap(u_a) + (1+nu)/2 * d_a(div u) - (1/ell^2) u_a = 0,
    1/ell^2 = k/(h Ebar).

Hence the only intrinsic length is

    ell = sqrt(h Ebar / k) = sqrt( h E / (k (1-nu^2)) ).      (decay length)

Far from the crack line the field decays over `ell`; this is the crack-induced
wake decay length. A uniform in-plane stretch `u_1 = gamma x_1` is excluded
because the foundation term `(k/2)u_1^2` would diverge over the infinite
`x_1` extent, so the bounded wake has `u_1` independent of `x_1`.

## 2. Far-ahead and far-behind states

### Case A (cracks on x2 = mH, tips at x1 = 0)
* Far ahead (x1 -> +inf): pristine prestressed state, `u = 0`.
* Far behind (x1 -> -inf): fully cracked wake, periodic in x2 with period H.
  Consider one strip `0 < x2 < H`; its edges `x2 = 0, H` are crack faces.
  |x2|-independent solution of the reduced system:
  - `u_1 = 0` (its ODE is `(1-nu)/2 u_1'' - u_1/ell^2 = 0` with `u_1'(0)=u_1'(H)=0`).
  - `u_2'' = u_2/ell^2`, giving
        u_2(x2) = B sinh((x2 - H/2)/ell),
        B = -sigma0 ell / (Ebar cosh(H/(2 ell))).
  - `sigma_22 = sigma0 + Ebar u_2' = sigma0[1 - cosh((x2-H/2)/ell)/cosh(H/(2ell))]`,
    `sigma_11 = sigma0 + Ebar nu u_2'`, `sigma_12 = 0`.
  - Strip-edge conditions (traction-free crack faces):
        `u_2' = -sigma0/Ebar`, `u_1' = 0`, equivalently `sigma_22 = sigma_12 = 0`
    at both `x2 = 0` and `x2 = H`.

### Case B (existing infinite cracks on x2 = 2mH for all x1; new semi-infinite
### cracks on x2 = (2m+1)H, x1 < 0, tips at x1 = 0)
* Far ahead: the equilibrium wake of the *old* cracks, i.e. the Case-A wake
  with period `2H` (strip half-width `H`, so `beta = H/ell`). Old cracks are
  traction-free for all `x1`; new tips have not yet passed.
* Far behind: the new cracks have split each old strip of width `2H` into two
  strips of width `H`, each in the Case-A wake state with half-width `H/2`
  (`beta = H/(2ell)`). Period is `2H`; there is exactly one advancing tip per
  period.

## 3. Modes permitted

Both geometries are symmetric about the (periodic) crack planes: `u_2` is odd
and `u_1` even about each crack plane. The applied load is equibiaxial normal
stress. The shear traction on the crack planes is identically zero, so the
crack-tip field is pure Mode I:

    K_II = 0,   G = K_I^2 / E.

Thus both cases permit only Mode I opening; `K_I` is reported below.

## 4. Energy release (per advancing tip)

`Delta Pi = int[ h sigma0 eps_gg + (h/2) Delta sigma_{ab} eps_{ab}
+ (k/2) u_a u_a ] dA` (PROBLEM.md). The loading is displacement-controlled
(`u -> 0` far ahead), so the release per unit advance is `-d(Delta Pi)/da`.
In steady state a strip of width `da` ahead (energy density `w_ahead = 0`,
since `u=0`) transitions to the wake state (`w_behind`). Hence, per period and
per unit advance,

    release = w_ahead - w_behind = -w_behind.

Wake energy per unit x1-length (strip width `W`, `beta_W = W/(2ell)`), obtained
by direct integration using `u_2'' = u_2/ell^2`:

    W_behind(W) = -(h sigma0^2 ell / Ebar) tanh(W/(2 ell)).

(Independent quadrature confirms this to 1e-16; see checks.) The linear term
`int h sigma0 eps_gg dx2` integrates to zero over a symmetric strip, and the
quadratic terms give the `tanh`.)

* Case A: one strip of width `H` per period, one tip per period:
      `G_A h = -W_behind(H)`, therefore
      G_A = (sigma0^2 ell / Ebar) tanh(H/(2 ell)).
* Case B: ahead the strip has width `2H`; behind it is two strips of width `H`:
      `G_B h = -W_behind(2H) + 2 W_behind(H)`, therefore
      G_B = (sigma0^2 ell / Ebar) [ 2 tanh(H/(2 ell)) - tanh(H/ell) ].

The thickness normalization `G = release/(h * da)` and the one-tip-per-period
count are exactly the ones requested ("release per unit in-plane advance per
tip divided by h").

## 5. Normalized results

With `beta = H/(2 ell)`:

    gamma_A = G_A Ebar/(sigma0^2 ell) = tanh(beta),
    gamma_B = G_B Ebar/(sigma0^2 ell) = 2 tanh(beta) - tanh(2 beta),
    kI_A = K_I,A/(sigma0 sqrt(ell)) = sqrt((1-nu^2) tanh(beta)),
    kI_B = K_I,B/(sigma0 sqrt(ell)) = sqrt((1-nu^2)(2 tanh(beta) - tanh(2 beta))).

Sampling at `ell/H = 0.5, 1, 1.5, 2` (i.e. `beta = 1, 0.5, 1/3, 0.25`),
`nu = 0.3`, `1-nu^2 = 0.91`:

| ell/H | gamma_A | gamma_B | kI_A | kI_B |
|:-----:|:-------:|:-------:|:----:|:----:|
| 0.5 | 0.7615941560 | 0.5591607318 | 0.8324966558 | 0.7133276007 |
| 1.0 | 0.4621171573 | 0.1626401586 | 0.6484802334 | 0.3847109880 |
| 1.5 | 0.3215127375 | 0.0602425297 | 0.5409034952 | 0.2341382114 |
| 2.0 | 0.2449186624 | 0.0277201675 | 0.4720974293 | 0.1588249114 |

(`GA, GB` follow from `G = gamma * sigma0^2 ell / Ebar`.)

## 6. Critical reference stress and limits

For toughness `Gc > 0`, advance starts when `G = Gc`:

    sigma0_c(A) = sqrt( Gc Ebar / (ell tanh(beta)) ),
    sigma0_c(B) = sqrt( Gc Ebar / (ell (2 tanh(beta) - tanh(2 beta))) ).

Limits at fixed `h, E, nu, k` (so `ell` fixed):

* `H/ell -> 0` (`beta -> 0`):
  - Case A: `G_A ~ sigma0^2 H/(2 Ebar)`, so `sigma0_c(A) ~ sqrt(2 Gc Ebar / H) -> inf`.
  - Case B: `G_B ~ sigma0^2 H^3/(4 ell^2 Ebar)`, so
    `sigma0_c(B) ~ sqrt(4 Gc ell^2 Ebar / H^3) -> inf` (faster).
  Meaning: with closely spaced cracks each tip relieves only a thin slice, so
  the available energy per tip is small; Case A grows like `H`, Case B like
  `H^3` because the new cracks only split an already relaxed old-crack wake.
* `H/ell -> inf` (`beta -> inf`): both `gamma_A, gamma_B -> 1`, hence
  `G_A = G_B -> sigma0^2 ell / Ebar` and `sigma0_c = sqrt(Gc Ebar/ell)`.
  Meaning: isolated (well-separated) cracks, each behaving as a single
  semi-infinite crack in a film on a foundation; the finite spacing and the old
  cracks become irrelevant.

## 7. Checks (executed, see check_output.txt; code in compute.py / checks.py)

1. Wake energy: direct numerical quadrature of the energy density vs the
   closed form `-(h sigma0^2 ell/Ebar) tanh(H/(2ell))`; max relative error
   1.5e-16 (tol 1e-9). Detects an algebraic slip in the wake energy.
2. Wake ODE residual `u_2'' - u_2/ell^2`; < 7e-7 (tol 1e-4). Detects an error
   in the decay length / ODE sign.
3. Traction-free strip edges `sigma_22 = 0`; exactly 0 (tol 1e-10). Detects a
   wrong boundary condition or normalization.
4. Independent numerical BVP solve (scipy `solve_bvp`) of the wake ODE vs the
   closed form; < 1e-10 (tol 1e-6). Detects a wrong general solution.
5. Integration-by-parts identity `int[(u_2')^2 + u_2^2/ell^2] = [u_2 u_2']`;
   rel 1e-16 (tol 1e-8). Detects an error in the energy reduction.
6. Small-`H` limits: `gamma_A ~ beta`, `gamma_B ~ 2 beta^3` (verified at
   `H=1e-3, 1e-2`). Detects a wrong coefficient in the thin-strip limit.
7. Large-`H` limits: `gamma_A, gamma_B -> 1`. Detects a wrong isolated-crack
   limit.
8. `G = K_I^2/E` consistency (`kI = sqrt(gamma (1-nu^2))`); exact (tol 1e-12).
   Detects an inconsistency between the energy and K normalizations.

Shared assumption: checks 1-8 all use the same 1D reduction/energy framework.
They are consistency checks, not an independent solution route. Result-set
`results.json` lists these with names/outputs/tolerances.

## 8. What is supported / what needs more

Supported by this reduced model:
- The wake decays over `ell = sqrt(h Ebar/k)`.
- Pure Mode I; closed-form `G_A(H)`, `G_B(H)` and their thin/thick limits.
- The ordering `G_A > G_B` for all finite `H` (since `2 tanh beta - tanh 2beta
  < tanh beta` for `beta > 0`): pre-existing cracks reduce the energy available
  to the new tips, and the effect is strongest at small `H`.
- `sigma0_c(B) > sigma0_c(A)`: the pre-cracked geometry needs a higher far
  stress to advance the new tips.

Requires additional information / validation (unresolved):
- An independent 2D solution (e.g. a J-integral on a boundary-element or
  finite-difference solve) was not completed within the budget; the energy
  coefficient therefore rests on the 1D reduction (see PROVISIONAL_REPORT.md
  and results.json "unresolved"). An error there would change the numerical
  prefactors but not the qualitative limits/ordering.
- The model assumes bonding is maintained, plane stress, linear elasticity and
  a fixed reference prestress during advance; deviations (interfacial
  debonding, large deflection, foundation yielding) are outside it.

## Reproduction

Single command (see REPRODUCE.md):

    timeout 180 python3 compute.py && timeout 180 python3 checks.py

Outputs: `results_analytic.json`, and the checks table printed to stdout
(saved to `check_output.txt`). The final machine-readable answer is
`results.json`.
