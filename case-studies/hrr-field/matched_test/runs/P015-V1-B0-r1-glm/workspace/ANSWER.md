# Answer: slip and strain contours near a stationary slit (mode III, power law)

## Requested numbers (mm, MPa)

| Quantity | Value |
|---|---|
| 1. Δw(0.0200 mm) | **+1.3181 × 10⁻³ mm** |
| 1. Δw(0.0800 mm) | **+1.7937 × 10⁻³ mm** (ratio (0.08/0.02)^(2/9) = 1.3606 ✓) |
| 2. τ_rz(0.0200 mm, 60°) | **+38.59 MPa** |
| 2. τ_θz(0.0200 mm, 60°) | **+138.0 MPa** |
| 3. r where \|γ\| = 0.0100, θ = 0° | **4.689 × 10⁻² mm** |
| 3. θ = 90° | **2.507 × 10⁻² mm** |
| 3. θ = 180° (upper face, limit) | **1.340 × 10⁻² mm** |

Signs: Δw > 0 (upper face moves +z as specified); both stress components positive on the
positive basis (e_r, e_θ), as in the elastic mode III limit (sin(θ/2), cos(θ/2) at 60°).

## Model chosen and why

Anti-plane (mode III) quasistatic equilibrium with the stated pure power law reduces to
div τ = 0 with τ = C|γ|^(1/n−1)γ, γ = ∇w, C = τ0γ0^(−1/n) (a "p-Laplacian" with p = 1/n).
The leading near-tip field is separable, w = A r^α f(θ). Because a strain energy exists,
loading is proportional and monotonic (fixed pattern), and the flanks are traction-free and
parallel to x, J is path-independent; for a separable crack field J ∝ r^(1−p(n+1)), so a
finite, nonzero J (**the given J = 0.120 N/mm**) forces the stress exponent
**p = 1/(n+1) = 2/9**, i.e. α = 1/(n+1) = 2/9, w ∝ r^(2/9), γ ∝ r^(−7/9), τ ∝ r^(−2/9)
— the mode III analogue of the HRR field. The angular function solves the nonlinear ODE

  (h^(k−1) f′)′ + (n/(n+1)) α f h^(k−1) = 0,  h = sqrt(α²f² + f′²),  k = 1/n,

with f(0) = 0 (the imposed antisymmetry w(r,−θ) = −w(r,θ)) and normalization f′(0) = 1;
traction-free faces need τ_θ(r,±π) = 0 ⇔ f′(π) = 0. I integrated the ODE (scipy solve_ivp,
rtol 1e-12) and fixed the single amplitude from J: with W = (n/(n+1))τγ for the pure power
law on the proportional path, J = C A^((n+1)/n) I_θ, I_θ = 2.01960, giving
A = 9.25609 × 10⁻⁴ mm^(7/9). That is how the loading fixes the magnitude: J is the sole
carrier of remote geometry/load into the leading field (all requested points have
r ≤ 0.08 ≪ 0.30 mm, where the problem states the leading field dominates).

## Assumptions changed

None of the stated assumptions were changed. Two additions: (i) normalization f′(0) = 1
(scale absorbed into A, then fixed by J — invariant, since scaling w by λ rescales
J by λ^((n+1)/n)); (ii) deformation-theory energy W and path-independent J were used, which
is exact here because inside the separable zone a fixed material point sees stress growing
proportionally (fixed direction, growing amplitude).

## Verification (`solve.py`, `eigen_check.py`)

1. **Boundary residual**: |τ_θz(r,±π)|/|τ| = |f′(π)|/(αf(π)) = 6.4 × 10⁻¹⁶ (faces
   traction-free, not imposed but satisfied — and perturbing the exponent by 10⁻³ breaks it,
   f′(π) ≈ −1.7 × 10⁻³, confirming the eigenvalue).
2. **Energy flux**: J recomputed from the field on closed rectangles (half-widths 0.05 and
   0.12 mm) with Cartesian ∂w/∂x by finite differences and proper polar-basis handling:
   J = 0.120000 N/mm on both — matches the supplied flux and shows path independence.
3. **Equilibrium**: div τ = 0 checked in Cartesian components by finite differences at
   several points off the faces: relative residual ≤ 1 × 10⁻⁹; ∂w/∂x FD vs analytic γ_x
   agrees to 1 × 10⁻¹¹.
4. **Known limiting solution**: with n = 1 the solver returns f = 2 sin(θ/2) (err 3.7 × 10⁻¹³),
   I_θ = π exactly, and the amplitude relation reproduces the Westergaard mode III result
   K_III = √(2GJ) to 9 × 10⁻¹⁴ (G = τ0/γ0).

Note: the provided `knowledge/checks/crack_tip_checks.py` targets in-plane (mode I/II, J2
plane strain/stress) fields and cannot ingest anti-plane shear, so I ran the analogous
checks myself (faces, equilibrium, compatibility via ∂w/∂x, J consistency, elastic limit);
the topic's p = 1/(n+1) result for cracks was reproduced independently.

## Missing / uncertain

- Small strain is marginal at the reported radii (|γ| up to ~2% at r = 0.02 mm) and fails as
  r → 0 (γ ∝ r^(−7/9)); the field is an intermediate asymptotic (inner limit set by
  continuum/small-strain breakdown, outer by r ≈ 0.30 mm per the problem).
- Higher-order (non-singular) terms were neglected per the problem statement; they set the
  true error of the leading-term values at the larger radii.
- n = 3.5 gives a smooth field (no sectoring), so the separable form is appropriate; nothing
  here addresses blunting or face contact (faces open: Δw > 0).

**Do these results alone determine the load at which the slit starts to advance?**
No — they characterize the near-tip field and its energy flux J, but advance requires a
fracture criterion (e.g., a critical J, critical crack-tip slip, or a process-zone/cohesive
strength condition), which is a separate material property not supplied here.
