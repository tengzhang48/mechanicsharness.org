# Answer — slip and strain contours near a stationary slit (mode III, power-law)

## The seven requested numbers

| # | Quantity | Value |
|---|---|---|
| 1 | Δw(r = 0.0200 mm) | **1.3181 × 10⁻³ mm** (≈ 1.32 µm) |
| 1 | Δw(r = 0.0800 mm) | **1.7937 × 10⁻³ mm** (≈ 1.79 µm) |
| 2 | τ_rz(r = 0.0200 mm, θ = 60°) | **+38.59 MPa** |
| 2 | τ_θz(r = 0.0200 mm, θ = 60°) | **+138.01 MPa** |
| 3 | r where \|γ\| = 0.0100 on θ = 0° | **0.04689 mm** |
| 3 | r where \|γ\| = 0.0100 on θ = 90° | **0.02507 mm** |
| 3 | r where \|γ\| = 0.0100 on θ = 180° (face limit) | **0.01340 mm** |

All requested positions lie at r ≤ 0.08 mm ≪ 0.30 mm, inside the stated validity of the leading field.
Numerical uncertainty of the tabulation/quadrature ≲ 2 × 10⁻⁷ relative (grid-refinement check).

## Governing model and why

- **Problem type.** Anti-plane (mode III) shear of a homogeneous body, w(x,y) only, γ = ∇w,
  equilibrium div τ = 0, faces traction-free (τ_θz = 0 at θ = ±π), antisymmetry w(r,−θ) = −w(r,θ).
  The specified law γ = (γ₀/τ₀)(τ/τ₀)ⁿ⁻¹ τ inverts to τ = τ₀(γ/γ₀)^{1/n} γ̂, so w solves a
  scalar nonlinear Laplace (p-harmonic) equation with p = (n+1)/n = 9/7. Loading is quasistatic,
  monotonic, fixed pattern ⇒ an energy density exists and J is path-independent (faces
  traction-free and parallel to x, no body forces — the topic's conditions).
- **Leading field.** A separable near-tip field w = r^α_w F(θ), τ ~ r^{−p} carries a finite,
  nonzero J only for 1 − p(n+1) = 0, i.e. **p = α_w = 1/(n+1) = 2/9** — the mode III analogue of
  the HRR exponent (mode I/II topic result applied to mode III). The angular functions follow
  uniquely (up to one amplitude) from the equilibrium ODE + face + symmetry conditions; they have
  a closed form obtained through the stress-function angle ψ, t = tan ψ ∈ [0, ∞):
  - θ(t) = atan t + atan(n t)  (θ: 0 → π as t: 0 → ∞; this map satisfies θ(π)=0-condition exactly)
  - T(t) = a [(1+t²)/(1+n²t²)]^{1/(2(n+1))}, a = n/(n+1), with T(0) = a, T(π) = a n^{−1/(n+1)}
  - τ_rz = A r^{−1/(n+1)} T sinψ,  τ_θz = A r^{−1/(n+1)} T cosψ
  - \|γ\| = c Aⁿ r^{−n/(n+1)} Tⁿ,  c = γ₀/τ₀ⁿ,  γ ∥ τ
  - w = D r^{1/(n+1)} F(θ),  D = c Aⁿ,  F(θ) = ∫₀^θ Tⁿ cos ψ′ dθ′, F(0)=0 (no rigid translation;
    antisymmetry gives w(r,0)=0), F(π) ≈ 0.70476
- **How the loading fixes the magnitude.** J scales as A^{n+1}: evaluating the given
  J = ∫(W N_x − t_z w,ₓ) ds (Cartesian components; W = (n/(n+1)) τγ for the pure power law) on any
  circle in the leading-field zone gives J(A=1) = 2.8466 × 10⁻¹⁰ (A in MPa·mm^{2/9}), hence
  **A = (0.120/2.8466×10⁻¹⁰)^{1/4.5} = 82.534 MPa·mm^{2/9}**, D = cAⁿ = 2.2307 × 10⁻³ mm^{7/9}.
  For reference: Δw(r) = 2D r^{2/9} F(π); r_{γ}(θ) = (D T(θ)ⁿ / 0.01)^{9/7}.

## Assumptions I changed (and why)

None of substance — the statement was used as given (pure power law with no elastic term,
deformation-type law justified by the stated monotonic proportional loading; J = 0.120 N/mm taken
as fully carried by the leading mode III field, which the statement grants for r < 0.30 mm).
Two judgement calls, both favourable: (i) small strain is used although \|γ\| → ∞ as r → 0 — at
the requested positions the engineering shear is only 1.0–1.5%, comfortably small; (ii) face
opening/contact is not an issue for antisymmetric mode III loading (faces slide, no interpenetration).

## Verification (all in `solve.py`, output below)

1. **Energy flux**: J computed from the field on circles R = 0.005, 0.02, 0.08, 0.20, 0.30 mm
   gives 0.12000000 N/mm at every radius (r-independence follows from the exponent; agreement of
   the absolute value with the target fixes the amplitude).
2. **Independent ODE integration**: solving the angular equilibrium BVP
   (T^{n−1}g′)′ + [n/(n+1)²] T^{n−1} g = 0, g(0)=−1, g′(0)=0 by stiff IVP integration reproduces
   the closed form to 1.6 × 10⁻¹⁵ and returns g(π) = 1.6 × 10⁻¹⁶ — confirming the eigen-consistency
   (traction-free faces) that the closed form encodes exactly. (This check caught and fixed a real
   algebra bug in my first closed-form T.)
3. **n = 1 regression** (law becomes LEFM, G = τ₀/γ₀): angular functions reduce to
   sin(θ/2), ½sin(θ/2), ½cos(θ/2) to ≤ 5 × 10⁻⁷, and the same J code yields K = A√(π/2),
   J = K²/(2G) exactly (I_J = π/4 to 3 × 10⁻¹³) — i.e. w = (2K/G)√(r/2π) sin(θ/2).
4. **PDE/BC residuals of the full field**: Cartesian finite-difference equilibrium residual
   ≤ 1.7 × 10⁻⁷ of scale; τ_θz on the faces ≤ 3 × 10⁻⁶ MPa (vs ~10² MPa field) — zero in closed
   form, the residual is tabulation truncation; constitutive check τ = τ₀(γ/γ₀)^{1/n} exact to
   printed digits; ∇w by finite differences matches the analytic γ exactly; the reported contour
   radii back-substitute to \|γ\| = 0.010000.

## Missing / uncertain

- The amplitude rests on J-dominance being valid to 0.30 mm as stated; no higher-order terms
  (needed only if the outer geometry mattered) were included.
- Strains grow without bound as r → 0, so the innermost reported point (0.0134 mm, 1% shear) is
  still fine, but points much closer to the tip would leave small strain.
- **Do these results determine the load at which the slit starts to advance?** No — they fix the
  near-tip fields at the current load only, and onset of growth requires a separate material
  failure criterion (e.g. a critical J or a critical slip at a characteristic distance), which is
  not contained in the data given.
