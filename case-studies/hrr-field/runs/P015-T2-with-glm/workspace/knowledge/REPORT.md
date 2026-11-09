# REPORT.md

## What was reproduced

- **Singularity exponent** `s=(2n+1)/(n+1)`. Independently recovered from the
  angular eigenvalue problem (both plane stress eq. 13 and plane strain eq. 31)
  by shooting on `(s, f''(0))`; agreement to `<10^{-6}` for `n=3,5,9,13`. The
  angular ODEs were integrated with `scipy.integrate.solve_ivp`; the highest
  derivative `f''''` is obtained in closed form (derived once with sympy, then
  implemented with numpy only).
- **Table 1 `I`** for `n=3,5,9,13`, plane stress and plane strain, via the paper's
  eq. (24): all eight within **0.3%** (ours 3.855/3.408/3.030/2.871 and
  5.507/5.024/4.602/4.401).
- **Table 1 `E_n`** via eq. (34): all eight within **0.9%**.
- **Table 1 `(π/I)^{1/(n+1)}`**: `<3×10^{-3}`.
- **Qualitative angular stress distributions** for plane stress and plane strain
  match Figs. 3a/b and 5a/b (checked by eye against the page images).
- **Linear-elastic special solution**: the first-principles circular J integral
  `elastic_J_numeric` on the mode-I biharmonic field (eq. 40) returns
  `J=K_I²/E` to quadrature accuracy and agrees with the closed form
  `elastic_J_plane_stress` (test `test_elastic_J_plane_stress_export_matches_numeric`).
  The integrand includes the polar-basis (angular-derivative) terms; an earlier
  version omitted them and used the wrong coefficient, and its agreement was a
  cancellation, not a validation.

## What was not reproduced and why

- **Piecewise-linear results (eqs. 40–45)** and the resulting amplitude factors
  `sqrt(E_tan/E)` (plane stress, eq. 42: `K=(1/√2)√(E_tan/E)σ^∞`, so the stresses
  are reduced as `E_tan` falls) and the plane-strain expression (eq. 45) were
  read but not turned into code: they require the yielded-region constitutive
  (eq. 35) and a separate J evaluation. `taken` from the paper. (The earlier
  `sqrt(E/E_tan)` was the reciprocal of what eq. (42) prints.)
- **The plane-strain/plane-stress tensile-stress ratio (eq. 50, Fig. 7)** was not
  computed numerically. It is a direct algebraic combination of quantities we did
  reproduce (`I`, `σ̃_θ(0)`), so it could be added, but was left `taken`.
- **Absolute amplitudes `K`** are `taken` (they follow from the reproduced `I` by
  J-matching); no independent far-field matching was attempted (the paper itself
  bypasses it).

## Resolved inconsistencies / remaining observations

1. **Eq. (22), traction-work line — printed prefactor error (resolved).** The image
   prints `K^{n+1} r^{(n-1)(s-2)}`, omitting `α` and writing `n-1`; dimensional
   consistency with the `W` line of the same equation and with (23)–(24) requires
   `αK^{n+1} r^{(n+1)(s-2)}`. This does not affect Table 1, since the printed eq.
   (24) I implements the correct scaling. See `EXTRACT.md`.
2. **The former "first-principles J discrepancy" was an artifact (resolved).** The
   earlier `I≈5.65` (plane stress `n=3`) came from dropping the derivatives
   `ũ_r'`, `ũ_θ'` of the displacement angular functions; a correct Cartesian route
   gives `I=3.8549663167`, matching Table 1. Eq. (24) as printed is correct.
3. **`E_n(n=9, plane strain)`**: image reads `1.48`, `paper.txt` gives `1.43`.
   Our computation (`1.485`) supports `1.48`.
4. **Cause of the Table 1 gaps is unestablished.** Agreement is only to
   three/four significant figures, but refining the present quadrature/ODE changes
   `E_n` by ~1e-8, far below the worst gap (0.85%, plane strain `n=5`). The
   difference from the printed table is therefore *not* an integration error of the
   present code, and its historical cause is not established here.

## Extraction items I am least sure of

- The exact typography of the displacement line in eq. (24): `.` (dot) vs nothing
  on `ũ_θ`. I interpreted the dotted quantities as `d/dθ` and used the
  strain–displacement relations to supply `ũ_θ` (`3σ̃_e^{n-1}σ̃_{rθ}-ũ_r')/m`); this
  choice is what reproduces Table 1.
- The plane-strain `E_n` row digits (`3` vs `8` glyphs); resolved by matching the
  computed values.
- Eq. (42)/(45) amplitudes (piecewise-linear) — only partially legible and not
  verified.

## Verification summary

`python3 -m pytest test_hrr.py` → **70 tests** (64 in round 2, 5 added in
round 3, 1 added in round 4): Table 1 checks (I, E_n, amplitude column), eigenvalue exponent checks,
tabulated energy/shape identity, an independent Cartesian finite-difference J
route for `I`, the elastic special solution (self-contained first-principles J
integral plus the `elastic_J_plane_stress`/`elastic_J_numeric` exports, with the
**returned** pointwise integrand checked against full Cartesian displacement
gradients), coverage of every exported function (`effective_stress_coeff`,
`amplitude_K`, `elastic_angular_stress_function`, `elastic_J_integrand`,
`elastic_displacement_cartesian`, `elastic_J_plane_stress`, `elastic_J_numeric`),
the J-scaling exponent identity, the corrected `amplitude_K` documentation and
plane-strain J reconstruction, `I`'s field-scale (`n+1`) normalisation, the
`1/E` scaling of the Cartesian displacement helper, and two planted-error checks
(wrong exponent; wrong effective stress). All pass.

Round-3 additions (see `CORRECTIONS.md`): the elastic component test now reads
the actual returned integrand (`resid = max|integrand - (W cos - traction)| < 1e-8`),
so restoring the old derivative expression in the returned value fails it
(`resid = 1.4`); `test_elastic_displacement_cartesian_modulus_scaling` calls the
joint at `E != 1`; `test_amplitude_K_documented_reconstruction` calls
`amplitude_K` exactly as documented and pins
`K = (J_dimless/(alpha I))^(1/(n+1))` for both constraints;
`test_integral_I_scales_as_field_to_n_plus_1` binds the stated normalisation
scope.

Round-4 addition (see `CORRECTIONS.md`): `test_documented_signatures_and_usage_statements`
calls each of the 12 documented signatures exactly as written, checks the
`sol`/`I` argument-binding claims against the real API, and reads `ENTRY.md`, so
it fails on the uncorrected prose (G1 of the round-3 verification) and passes on
the corrected text.
