# Slip and strain contours near a stationary slit (mode III, power-law)

## Requested numbers (mm, MPa; J = 0.120 N/mm, n = 3.5)

| # | Quantity | Value |
|---|----------|-------|
| 1a | Δw(r = 0.0200 mm) | **1.318 × 10⁻³ mm** (0.0013181) |
| 1b | Δw(r = 0.0800 mm) | **1.794 × 10⁻³ mm** (0.0017937) |
| 2a | τ_rz(r = 0.0200 mm, θ = 60°) | **+38.59 MPa** |
| 2b | τ_θz(r = 0.0200 mm, θ = 60°) | **+138.0 MPa** |
| 3a | r where \|γ\| = 0.0100 at θ = 0° | **0.0469 mm** (0.046894) |
| 3b | r where \|γ\| = 0.0100 at θ = 90° | **0.0251 mm** (0.025066) |
| 3c | r where \|γ\| = 0.0100 at θ = 180° (face limit) | **0.0134 mm** (0.013398) |

Δw > 0 (upper face moves +z, consistent with the stated loading); ratio Δw(0.08)/Δw(0.02) = 4^{1/(n+1)} = 1.3617 exactly.

## Governing model and how the loading fixes its magnitude

Anti-plane shear (u_z = w only), equilibrium div τ = 0, compatibility curl γ = 0 with
γ = ∇w and the pure power law γ = (γ₀/τ₀ⁿ) τ^{n−1} τ. With a stress function
ψ = B r^{1−s} Φ(θ) (τ_r = ψ_θ/r = B r^{−s}Φ′, τ_θz = −ψ_r = −(1−s) B r^{−s}Φ),
separability gives the nonlinear eigenproblem

  d/dθ[ T^{n−1} Φ′ ] + (1−sn)(1−s) Φ T^{n−1} = 0,  T² = Φ′² + (1−s)²Φ²,
  Φ(0) = −1 (sign so upper face moves +z), Φ′(0) = 0 (τ_rz odd, w odd in y),
  Φ(π) = 0 (traction-free faces ⇒ τ_θz = σ_yz = 0 on θ = ±π).

The exponent s = 1/(n+1) (= 2/9) is required for the field to carry finite,
radius-independent energy flux (J-dominance: W ~ τ^{n+1} ~ r^{−s(n+1)}); the numerical
shooting eigenvalue reproduces s* = 0.2222222222 with face residual Φ(π) = 1.5×10⁻¹⁶,
confirming it. This is the mode-III HRR field of a pure power-law material.

The loading magnitude enters only through J. Evaluating J on a circle for this field,
J = (γ₀/τ₀ⁿ) B^{n+1} I with I = ∫T^{n−1}[(n/(n+1))T²cosθ − Φ′²cosθ − (1−s)ΦΦ′sinθ]dθ
= 0.651800 (n = 3.5), giving the amplitude B = (J τ₀ⁿ/(γ₀ I))^{1/(n+1)} = 82.534 MPa·mm^{2/9}.
Then: τ_rz = B r^{−s}Φ′, τ_θz = −(1−s)B r^{−s}Φ, τ = B r^{−s}T, |γ| = γ₀(τ/τ₀)ⁿ,
w = (γ₀/τ₀ⁿ)(Bⁿ/s) r^s T^{n−1}Φ′ (no rigid translation), so
Δw = 2(γ₀/τ₀ⁿ)(Bⁿ/s) r^s [T^{n−1}Φ′]_{θ=π}, and |γ| = 0.01 at
r = [B T(θ)/(τ₀ (0.01/γ₀)^{1/n})]^{(n+1)}. T(0°) = 0.77778, T(90°) = 0.67671, T(180°) = 0.58878.

## Assumptions changed

None of the stated assumptions were changed. The only interpretive choices: J > 0 drives
the field forward, and the sign of Φ was fixed so the upper face moves in +z (as specified).

## Verification (quantitative)

1. **Face boundary residual:** τ_θz(r, ±π) = −(1−s)B r^{−s}Φ(π); computed Φ(π) = 1.5×10⁻¹⁶ ⇒ face traction ~10⁻¹⁴ MPa at r = 0.02 (zero to machine precision).
2. **Energy flux:** J recomputed by direct contour integration ∫(W cosθ − τ_rz γ_xz) r dθ of the reconstructed field gives 0.120000000 N/mm at r = 0.020, 0.080 and 0.200 mm (radius-independent, matches the supplied value).
3. **Known limiting solution:** setting n = 1 (μ = τ₀/γ₀ = 40 000 MPa) the whole pipeline reproduces the exact LEFM K-field: Δw = 4√(Jr/πμ), τ_rz = K sin(θ/2)/√(2πr), τ_θz = K cos(θ/2)/√(2πr), γ-contour radius K²/(2π(0.01μ)²), all to ≤ 3×10⁻¹⁶ relative error.

Code: `solve.py` (solve_ivp DOP853, rtol 1e-12; scipy quad for J).

## What is missing or uncertain

- The largest uncertainty is the truncation to the leading (separable) field, assumed negligible for r ≪ 0.30 mm per the problem statement; at the largest requested radius (0.08 mm ≈ 0.27×0.30 mm) higher-order terms from the far geometry could enter at the few-percent level (within the 2% tolerance was not independently checked, as no far-field geometry is given).
- ODE/quadrature tolerances contribute errors far below 10⁻⁸.
- **Do these results determine the load at which the slit starts to advance?** No — the fields are self-similar and singular at the tip (nonzero J with no threshold), so without an additional growth criterion (e.g., a critical J, CTOD, or cohesive strength) and the load–J relation of the specific specimen, the onset of advance cannot be determined.
