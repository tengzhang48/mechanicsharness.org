# ANSWER.md — Separation between unequal strips (moment-loaded, bridged)

## Results

| Quantity | Value |
|---|---|
| 1. Moment on each strip, M* | **0.2126 N·m** (0.212626) |
| 2. Advance from starter-cut root, L* | **10.80 mm** (10.8001) |
| 3. Opening at mid-zone, δ(L*/2) | **38.2 µm** (38.219) |

At this instant the trailing-edge connector (the oldest newly retained one, at
x = L*) first reaches δc and becomes inactive; the applied moment is at its
quasistatic maximum (peak of the R-curve), J* = 21.40 J/m².

## Method

**Bending stiffnesses** (plane strips, lateral contraction free, E as quoted):
D1 = E1 b h1³/12 = 0.04950 N·m², D2 = 0.288563 N·m²,
C_M = 1/D1 + 1/D2 = 23.6675 (N·m²)⁻¹.

**Mechanical description.** With no end forces, each arm between the bridged
zone and its free end carries constant moment M and zero shear (pure bending).
Far ahead the bonded part is long and straight; treating the intrinsic front
process as spatially negligible, the arms are clamped where they meet the
straight bonded part: w_i(0) = w_i'(0) = 0 at the front (x = 0). Let x run from
the front toward the starter-cut root (x = L). The newly separated faces
(0 < x < L) carry the closing connector traction t(δ) = σb(1 − δ/δc); the
starter-cut faces are traction-free, so at x = L the arms' moment is M and shear
zero: δ''(L) = C_M M, δ'''(L) = 0.

**Governing BVP.** The two Euler–Bernoulli arms obey D_i w_i'''' = −b t(δ),
δ = w1 + w2. Until the first connector goes inactive the law is linear, giving

    δ'''' = μ⁴ (δ − δc),   μ⁴ = b C_M σb / δc   (μ⁻¹ = 6.412 mm),
    δ(0) = 0,  δ'(0) = 0,  δ''(L) = C_M M,  δ'''(L) = 0.

(The softening law acts as a *negative*-stiffness foundation, hence exponential
rather than oscillatory modes.) Solved exactly with the basis
{cosh μx, sinh μx, cos μx, sin μx} (`solve.py`).

**Growth condition (Rice J-integral; Bao & Suo 1992 eqs (15), (22), (51), the
splitting-beam specimen).** Because the outer field is moment-only,

    J_far = M² C_M / (2b)   — independent of crack length,

and path independence of J across the bridged zone gives, at every instant of
quasistatic growth,

    M² C_M/(2b) = Γ0 + ∫₀^{δ_t} t dδ = Γ0 + σb δ_t − σb δ_t²/(2 δc),
    δ_t = δ(L).

δ(x) rises monotonically from the front, so the most-opened connector is the
trailing-edge one at x = L. **The event δ_t = δc** therefore saturates the
bridging work at ∫₀^{δc} t dδ = σb δc/2 = 20.0 J/m²:

    J* = Γ0 + σb δc/2 = 21.40 J/m²  →  M* = √(2 b J*/C_M) = 0.212626 N·m.

L* follows from the BVP with δ(L*; M*) = δc: **L\* = 10.8001 mm**, and the same
profile gives **δ(L\*/2) = 38.219 µm** (bending-only would give 50 µm; the
connectors hold the zone more closed).

**Path structure (why this is the quasistatic event).** Continuation from
initiation (M → √(2bΓ0/C_M) = 0.0544 N·m as L → 0) along
J_far = Γ0 + ∫t dδ shows M(L) rising monotonically and stably
(tangent stiffness G_M = C_M M/b − (∂δ_t/∂M) t(δ_t) > 0 throughout) to its
maximum exactly at δ_t = δc — structurally, dM/dL = 0 ⟺ t(δ_t) = 0. So the
first inactivation coincides with the peak moment; no connector loses activity
earlier (δ(x) ≤ δ(L) = δc inside the zone).

## Checks

1. **Energy consistency of the growth condition.** Computing the exact total
   potential energy Π (arms + connector springs, moment loading) of the BVP
   solutions, −∂Π/∂L = b[J_far − ∫₀^{δ_t} t dδ] holds to 3·10⁻⁶ (quadrature
   limit) at the event and elsewhere; at the event both equal bΓ0 = 0.035 J/m,
   and J_far = 21.4000 = Γ0 + σbδc/2 closes exactly.
2. **Independent ODE integration.** `scipy.integrate.solve_ivp` (adaptive RK,
   rtol 1e-12) with shooting on the tip curvature unknowns reproduces
   δ(L/2) = 38.2189 µm and δ(L) = 200.0000 µm — agreement with the closed-form
   solution to 6 digits.
3. **Limits.** L → 0: path M(L) → LEFM initiation value √(2bΓ0/C_M) = 0.0544
   N·m. σb → 0: BVP reduces to the clamped bending parabola δ = C_M M x²/2
   (connectors only reduce the opening, as they must).
