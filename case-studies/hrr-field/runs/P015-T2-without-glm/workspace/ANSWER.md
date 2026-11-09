# Answer

**1. p = 0.1562**  (σθθ(r,0) ∝ r⁻ᵖ; more precisely 0.15618, tolerance 0.0003)
**2. R = 0.3733**  (σrr(r,0)/σθθ(r,0); more precisely 0.37326, tolerance 0.002)

These describe the leading symmetric opening (mode-I-like) asymptote at the tip of the 120°-opening
notch (material wedge 240°, flanks at θ = ±2π/3) in the n = 4 power-law material.

## Model and boundary conditions (justification)

- **Plastic power law dominates at the tip.** The eigen-stress diverges like r⁻ᵖ, so σₑ → ∞ and the
  plastic strain ∝ (σₑ/σ₀)ⁿ⁻¹ overwhelms the elastic strain as r → 0. The leading asymptote is therefore
  the pure power-law (J₂ deformation-theory) material: E, ν, σ₀, α drop out — p and R depend only on
  n = 4 and the notch angle, and the remote load amplitude cancels (the eigenproblem is homogeneous in
  the stress amplitude; we normalize f(0) = 1 and choose the sign so σθθ(r,0) > 0).
- **Plane strain.** εzz = 0 with vanishing leading-order s_zz gives σzz = (σrr + σθθ)/2, and plastic
  incompressibility gives εrr = −εθθ at leading order.
- **Separable HRR-type ansatz.** σij = r⁻ᵖ σ̃ij(θ) from an Airy function Φ = r²⁻ᵖ f(θ):
  σ̃θθ = (2−p)(1−p) f, σ̃rr = (2−p) f + f″, σ̃rθ = −(1−p) f′, σ̃ₑ² = (3/4)(σ̃rr−σ̃θθ)² + 3σ̃rθ²;
  strains ∝ r⁻ⁿᵖ σ̃ₑⁿ⁻¹ s̃ij; displacements ∝ r¹⁻ⁿᵖ.
- **Governing ODE** (compatibility of those strains with u, v ∝ r¹⁻ⁿᵖ, with the correct 1/(1−np)
  factors from u = ∫εrr dr):
  (σ̃ₑⁿ⁻¹D)″ + np(2−np) σ̃ₑⁿ⁻¹D − 4(1−np)(σ̃ₑⁿ⁻¹T)′ = 0,  D = f″ + p(2−p)f,  T = −(1−p)f′.
- **BCs.** Traction-free flanks: f(±2π/3) = 0 and f′(±2π/3) = 0 (σθθ = σrθ = 0). Symmetry of the
  opening mode about θ = 0: f′(0) = f‴(0) = 0 (σrθ(0) = 0, uθ(0) = 0). Normalization f(0) = 1.
- **Solution.** Shooting from θ = 0 (DOP853, rtol 1e-12) with unknowns p and q = f″(0), imposing
  f(2π/3) = f′(2π/3) = 0. A full scan of p ∈ (0.01, 0.249) (p < 1/n = 0.25 required for decaying
  displacements) × q ∈ (−6, 2) finds a **unique** symmetric eigenmode:
  **p = 0.15617526, q = −1.2630790**, hence σ̃θθ(0) = (2−p)(1−p) = 1.55586,
  σ̃rr(0) = (2−p) + q = 0.58075, σ̃zz(0) = 1.06831, and
  **R = σ̃rr(0)/σ̃θθ(0) = 0.373262**.
- **Consistency of the field.** strains ∝ r⁻⁰·⁶²⁴⁷, displacements ∝ r⁰·³⁷⁵³ → 0 at the tip;
  energy density ∝ r⁻ᵖ⁽ⁿ⁺¹⁾ = r⁻⁰·⁷⁸¹ integrable (p(n+1) = 0.781 < 2, finite energy);
  min σ̃ₑ² = 0.71 > 0 on [0, 2π/3] (no sector degeneracy); p < 1/(n+1) = 0.2 as expected for a notch
  (blunting weakens the singularity relative to a crack), and R < 1 reflects relaxed triaxiality at the
  notch bisector (elastic 240° wedge gives R = 0.987; perfectly-plastic crack limit is 0.61).

## Quantitative verification

1. **Crack limit, exact values.** The same code with a 360° material wedge (crack, flanks at ±π) gives
   p = 0.250000000 (n = 3), **0.200000000 (n = 4)**, 0.100000000 (n = 9) — exactly 1/(n+1) (HRR),
   validating the ODE (with its 1−np factors), BCs and eigen-solver to 9 digits.
2. **Independent numerical method.** Global adaptive collocation two-point BVP solver
   (scipy `solve_bvp`, p as a free parameter; 9439 nodes, max relative residual 1.0e-11, boundary
   residual 5.7e-34): p = 0.1561752603 and R = 0.37326231 — differences from shooting
   Δp = 1.4e-11, ΔR = 2.4e-9. The crack case redone this way gives p = 0.2000000000.
3. **Numerics hygiene.** Shooting tolerance study rtol = 1e-9 … 1e-13: p, q, R identical to 10 digits
   (p = 0.1561752603, R = 0.37326231); boundary conditions met to ~1e-16; the trajectory satisfies the
   ODE pointwise (interior finite-difference check, rel ≈ 2e-7); equilibrium is satisfied identically
   by the Airy parametrization. (A global cosine-series collocation was also tried but abandoned: the
   eigenfunction has a steep 4th-derivative layer at the flank that a global basis cannot resolve.)

## Answers

| quantity | value |
|---|---|
| p in σθθ(r,0) ∝ r⁻ᵖ | **0.1562** (0.15618 ± 0.0003) |
| R = σrr(r,0)/σθθ(r,0) | **0.3733** (0.37326 ± 0.002) |

Independent of E, ν, σ₀, α and of the remote load level, as required.
