# ANSWER

## Results (mm–N–MPa units)

| Quantity | Value |
|---|---|
| **1. T** | **T = 4816 MPa** (t = T/E = 0.068807; 1% band 4768–4865) |
| 2. C_rr(60°) | **+3.942 × 10⁻³ mm^1/2** |
| 2. C_θθ(60°) | **+6.332 × 10⁻⁴ mm^1/2** |
| 2. C_γrθ(60°) | **+5.731 × 10⁻³ mm^1/2** (engineering shear = 2ε_rθ) |
| 3. C_y(90°) | **+9.139 × 10⁻³ mm^1/2** |
| 4. Radius interval on θ = 0 | **0.0100 mm ≤ r ≤ 0.0613 mm** (non-empty) |

Auxiliary: effective Poisson ratio ν* = 0.48624, effective stress-intensity factor K = 62.074 MPa·mm^1/2, K/√(2π) = 24.767 MPa.

## Model chosen and why

**Reduction of the constitutive law.** As r → 0, σ_e → ∞, so ⟨1 − σ_y/σ_e⟩ → 1 and the total-strain law becomes *exactly linear* in σ:

ε_ij = (1+ν)/E σ_ij − ν/E σ_kk δ_ij + (3/2) A s_ij,  A = 1/T − 1/E.

Collecting deviatoric and hydrostatic parts, this is an isotropic Hooke law with **effective moduli E\* = T** (since 1/E\* = 1/E + A = 1/T) and **ν\* = ν t + (1−t)/2 = ½[1 + (2ν−1)t]**, t = T/E. (Uniaxially: ε = σ/E + (σ−σ_y)A for σ > σ_y, i.e. a bilinear curve with post-yield tangent modulus T; the tip field sees only the tangent branch.)

**Field.** With stationary crack, proportional monotonic loading, and a total-strain law with potential W(σ) = σ:σ/(2E)-part + (A/2)⟨σ_e − σ_y⟩² (so ε = ∂W/∂σ), the leading near-tip equilibrium/kinematics problem is therefore the **plane-stress Williams mode-I field with (E\*, ν\*) = (T, ν\*)** — consistent with the reported √r scaling of ε_zz. ε_zz (free, with σ_zz = 0) comes out of the same reduction: on the tip region ε_zz = −[ν/E + A/2]⟨…⟩(σ_xx+σ_yy) → −(ν\*/T)(σ_xx+σ_yy), i.e. plain plane-stress Hooke with the effective constants; both the elastic and the plastic term feed ε_zz through the deviatoric part.

**Amplitude from J.** The corrections that distinguish the true law from the effective-linear law (terms ∝ σ_y/σ_e) are bounded in ε and O(r^−1/2) in W; multiplied by ds ~ r on a shrinking circular contour they contribute O(√r) → 0. Hence J equals its value on the effective-linear field, i.e. the Irwin relation for plane stress, J = K²/E\* = **K²/T**, so K = √(JT).

**Calibration of T.** On θ = 0: σ_xx = σ_yy = K/√(2πr), σ_xy = 0 ⇒ ε_zz = −(2ν\*/T)K/√(2πr) ⇒

C_zz = −2ν\*√(JT)/(√(2π) √T) = −0.005, with ν\* = ½[1 + (2ν−1)T/E].

Solving (brentq; LHS monotone in T, unique root) gives **T = 4816.47 MPa**, ν\* = 0.48624, K = 62.074 MPa·mm^1/2. Reproduction check: C_zz = −0.00500000.

**Parts 2–3.** Cartesian K-field stresses → strains via ε = [(1+ν\*)σ_ij − ν\*σ_kk δ_ij]/T, transformed to polar components at θ = 60° (equivalently σ_rr = C cos(u)(1+sin²u), σ_θθ = C cos³u, σ_rθ = C sin u cos²u, u = θ/2 — note σ_rθ(60°) ≠ 0). Displacements (no rigid translation): u = (K/(2μ\*))√(r/2π){cos(u)[κ\*−1+2sin²u], sin(u)[κ\*+1−2cos²u]}, μ\* = T/(2(1+ν\*)), κ\* = (3−ν\*)/(1+ν\*) = 1.69136; at θ = 90° the bracket is κ\*+1−1 = κ\*, giving C_y = K(1+ν\*)κ\* sin45°/(T√(2π)) = 9.139 × 10⁻³.

**Part 4 (screens on θ = 0, leading field).** On θ = 0 the stress is equibiaxial (σ, σ, 0) so σ_e = σ = K/√(2πr); principal strains are ε_xx = ε_yy = (1−ν\*)σ/T and ε_zz = −2ν\*σ/T, and since 2ν\* = 0.9725 > 1−ν\* = 0.5138, the largest magnitude is |ε_zz| = 0.005/√r (the measured C_zz itself).
- σ_e ≥ 10σ_y = 100 MPa ⇒ r ≤ K²/(2π·10⁴) = **0.061325 mm**;
- |ε_zz| ≤ 0.05 ⇒ 0.005/√r ≤ 0.05 ⇒ r ≥ **0.010000 mm** (ε_xx, ε_yy give the weaker bound r ≥ 0.00279);
- r ≤ 0.200 mm is inactive.
Interval: **[0.0100, 0.0613] mm**, non-empty.

## Assumptions changed

None of the stated assumptions were changed. I *added* the standard (and here forced) asymptotic reduction stated above — leading field = plane-stress Williams field at (T, ν\*), amplitude set by J = K²/T — which is a consequence of σ_e → ∞ at the tip, not a change of the given data. I also relied on the given statements that loading is proportional (⇒ deformation theory valid, J path-independent, W a potential) and that matching corrections are negligible for r ≤ 0.2 mm.

## Verification (beyond substituting into C_zz)

1. **J-integral on the leading field, done numerically** (independent of the C_zz measurement): I evaluated J = ∮(W N_x − t_i u_{i,x}) ds on circular arcs of radius 0.001, 0.01, 0.1 mm using the effective-linear W, tractions, and finite-difference displacement gradients: J = 0.80000000 at all three radii, matching the supplied J and K²/T to 8 digits. This validates both the field and the amplitude relation K = √(JT).
2. **Displacement field checked kinematically**: finite-difference strains from u(x,y) match the stress-derived strains to 7 significant digits at four sample points, and polar-derivative kinematics at θ = 60° (ε_rr = ∂u_r/∂r, ε_θθ = u_r/r + ∂u_θ/(r∂θ), γ_rθ = (1/r)∂u_r/∂θ + ∂u_θ/∂r − u_θ/r) reproduce 3.942048e-3, 6.331781e-4, 5.731131e-3 exactly — part 2 verified without using Hooke's law at all.
3. **Consistency of the reduction in the annulus**: evaluating the *full* nonlinear law on the leading stresses on θ = 0 changes ε_zz by 3.9% at r = 0.010 and 9.6% at r = 0.0613 (= σ_y/σ_e effect), confirming the σ_e ≥ 10σ_y screen bounds the constitutive-reduction error by ≈10% as intended.

## What the screens establish; what is missing

The screens delineate the annulus on θ = 0 in which the *leading field alone* is a self-consistent description: deeply yielded material (σ_e ≥ 10σ_y keeps the neglected σ_y/σ_e term ≤ 0.1, i.e. ≤ ~10% constitutive error), infinitesimal strain still defensible (principal strains ≤ 5%), and leading-field validity (matching corrections negligible for r ≤ 0.2 mm). They are defined consistency checks on the asymptotic field, **not** error estimates for the full finite-r solution (as the task states), and they say nothing about where the actual fields deviate in detail.

**Crack initiation is not determined.** The constitutive law plus J fix the stress/strain/displacement fields (and their coefficients), but no failure criterion (critical J or K, CTOD, critical strain/stress, damage model) is supplied; crack initiation and its location/mode therefore cannot be inferred from the given data.

## Uncertainty

- T inherits only the stated input precision (J, C_zz, E, ν given to ~3–4 digits); everything else is analytic, machine-precision.
- Within the required tolerances the answers are insensitive to the O(σ_y/σ_e) ≤ 10% constitutive correction (that correction affects the *true* field, not the requested leading-field coefficients, which are defined asymptotically).
