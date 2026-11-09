# Answer: slip and strain contours near a stationary slit (mode III, power law)

## 1. Requested numbers

**(1) Relative face slip Δw(r) = w(r,+π) − w(r,−π):**

| r (mm) | Δw (mm) |
|---|---|
| 0.0200 | **1.3181 × 10⁻³** |
| 0.0800 | **1.7937 × 10⁻³** |

**(2) Signed stresses at r = 0.0200 mm, θ = 60° (e_r = (cosθ, sinθ), e_θ = (−sinθ, cosθ)):**

- τ_rz = **+38.587 MPa**
- τ_θz = **+138.011 MPa**

(Both positive, as in the elastic mode III field ∝ sin(θ/2), cos(θ/2), for loading with the
upper face moving in +z.)

**(3) Radii where |γ| = 0.0100:**

| ray | r (mm) |
|---|---|
| θ = 0° | **0.046894** |
| θ = 90° | **0.025066** |
| θ = 180° (upper face, limit from inside) | **0.013398** |

All radii in (3) and all requested positions are < 0.30 mm, inside the stated validity region
of the leading field. Field quantities above use the eigenfunction normalization f(0) = −1
(results are normalization-independent, see §2).

## 2. Governing model and how the loading fixes its magnitude

**Model.** The only displacement is u_z = w and the loading is antisymmetric,
w(x,−y) = −w(x,y): this is pure **mode III (anti-plane shear)**. With monotonic quasistatic
loading, no unloading, and the given pure power law (no elastic term), the stress/strain state
at each point grows proportionally, so the stated total (deformation-theory) law
**γ = (γ₀/τ₀ⁿ) τⁿ⁻¹ τ** is the exact constitutive relation, a strain-energy density
W = (n/(n+1)) τ·γ exists, and Rice's J is path-independent (faces traction-free and parallel
to the x-axis, homogeneous material, no body forces). The near-tip field is therefore the
separable HRR-type field, written here through a stress function Φ = A r^λ f(θ):

- τ_rz = (1/r)Φ_,θ = A r^(−p) f′(θ),  τ_θz = −Φ_,r = −Aλ r^(−p) f(θ),  p = 1 − λ,
- |τ| = A r^(−p) T(θ),  T = √(f′² + λ²f²),
- |γ| = γ₀(τ/τ₀)ⁿ,  γ_i = (γ₀/τ₀ⁿ) τⁿ⁻¹ τ_i,
- w(r,θ) = (γ₀/τ₀ⁿ) Aⁿ (r^p/p) T(θ)ⁿ⁻¹ f′(θ)  (satisfies w_,r = γ_rz, (1/r)w_,θ = γ_θz; no rigid translation since w(0) = 0).

Compatibility (∇×γ = 0) gives one nonlinear eigenvalue ODE, with m = (n−1)/2, q = n(λ−1)+1:

  f″ = −[(n−1)λ² f f′² + λq (f′² + λ²f²) f] / [n f′² + λ² f²],

with f′(0) = 0 (antisymmetry ⇒ τ_rz = 0 ahead of the tip) and f(π) = 0 (traction-free face,
τ_θz = 0). For a crack, a finite nonzero J requires **p = 1/(n+1)**, i.e.
λ = n/(n+1) = 7/9 (topic §3); solving the ODE eigenvalue problem by shooting independently
returns λ = 0.777777777778, agreeing with n/(n+1) to 12 digits.

**Amplitude.** Substituting the field into the given J integral (on a circular arc the traction
is t_z = τ_rz and w_,x = γ_rz cosθ − γ_θz sinθ, so no polar-basis derivative pitfall arises for
the scalar w) gives

  J = (γ₀/τ₀ⁿ) Aⁿ⁺¹ Iₙ,  Iₙ = ∫₋π^π [(n/(n+1)) Tⁿ⁺¹ cosθ − Tⁿ⁻¹ f′(f′cosθ + λf sinθ)] dθ.

With J = 0.120 N/mm, n = 3.5, τ₀ = 80 MPa, γ₀ = 0.002:
Iₙ = 0.65179982 (normalization f(0) = −1; Iₙ and A rescale together, the field is invariant),
c ≡ γ₀/τ₀ⁿ = 4.36732 × 10⁻¹⁰ MPa⁻³·⁵, hence
**A = [J/(c Iₙ)]^(1/(n+1)) = 82.5342 MPa·mm^(2/9)**. Thus the single scalar J fixes the one free
constant of the leading field.

Evaluations: Δw(r) = 2w(r,π) = 2c Aⁿ (r^p/p) f′(π)ⁿ with f′(π) = 0.588778;
τ_rz, τ_θz from the formulas above; the |γ| = 0.01 contour, with
τ* = τ₀(0.01/γ₀)^(1/n) = 126.706 MPa, is r(θ) = [A T(θ)/τ*]^(n+1), using T(0) = λ = 0.777778,
T(π/2) = 0.676711, T(π) = f′(π) = 0.588778.

## 3. Assumptions changed

None. Every statement in the task was used as given: the specified pure power law (no elastic
term added), small strains, traction-free faces, the leading separable field at all requested
positions (the task states outer-geometry corrections are negligible outside r = 0.30 mm, and
all requested radii are ≤ 0.08 mm), no rigid translation (w(0) = 0), and the loading direction
only fixes the overall sign (chosen so the upper face moves in +z, giving Δw > 0).

## 4. Verification (quantitative)

All checks are in `solve.py` (output in `solution_output.txt`):

1. **Known limiting solution (n = 1).** With n = 1 the law is linear elasticity with
   μ = τ₀/γ₀. The same code returns λ = 0.5 exactly, f(θ) = −cos(θ/2) and f′ = ½ sin(θ/2)
   to ≤ 7 × 10⁻¹⁴, and I₁ = 0.785398163397 vs π/4 = 0.785398163398 — the exact elastic
   mode III crack field, with the J relation reducing to J = K²/(2μ) (checked algebraically:
   A = K√(2/π), cA²I₁ = K²/(2μ)). This validates the ODE, the eigenvalue search, the Iₙ
   quadrature and the amplitude relation end to end.
2. **Energy-flux calculation.** J was recomputed from the constructed n = 3.5 field on
   circular contours r = 0.001, 0.02, 0.08, 0.25 mm, using Cartesian components
   (W cosθ − τ_rz γ_xz) r dθ: all give 0.120 N/mm with relative error ≤ 5 × 10⁻¹⁴,
   confirming both the exponent (path independence) and the amplitude A.
3. **Boundary residual.** Face traction τ_θz(r, ±π) ≈ 3 × 10⁻¹⁴ MPa (exactly zero by the
   eigenfunction BC f(π) = 0, residual = −1.9 × 10⁻¹⁶ in f).
4. **Compatibility / displacement consistency.** (1/r)w_,θ vs γ_θz on 4001 points:
   max relative difference 2 × 10⁻⁷; w_,r vs γ_rz agrees to all printed digits.
5. **Equilibrium** τ_xz,x + τ_yz,y = 0 by central differences at four off-axis points:
   relative residual ≤ 1 × 10⁻⁷ of the local stress gradient scale.
6. **Independent recomputation (`verify2.py`).** The |γ| = 0.01 radii were re-obtained by
   root-finding |γ|(r) = 0.01 directly on the field (not by the closed-form inversion), and
   Δw was re-obtained by numerical quadrature of γ_rz along the face ray: agreement to
   ≤ 6 × 10⁻¹⁴. Doubling the angular grid (1e5 → 4e5 points) changes all seven reported
   numbers by ≤ 5 × 10⁻¹¹ relative.

## 5. Missing or uncertain

- The numbers rest on the leading separable term alone; higher-order (outer-geometry) terms
  are stated to be negligible for r ≤ 0.30 mm but were not quantified, so the true field may
  deviate near the outer requested radii by an unknown (small) amount.
- Small-strain theory is self-inconsistent sufficiently close to the tip (strains ≳ 0.01
  already at r ≈ 0.02 mm; a blunting/finite-strain zone of order J/τ₀ ≈ 1.5 × 10⁻³ mm is not
  modelled). The task instructs use of the leading small-strain field, so this does not change
  the answers but bounds their physical interpretation.
- Iₙ is quoted in the normalization f(0) = −1; any rescaling of f leaves A·T (and hence all
  reported fields) unchanged.
- The provided `checks/crack_tip_checks.py` targets in-plane (mode I/II) fields
  (S_rr, S_tt, S_rt), so it does not apply to this anti-plane problem; the equivalent
  mode-specific checks above were used instead.

## 6. Does this determine the load at which the slit starts to advance?

**No.** These results describe the stationary near-tip field at one given applied J; deciding
when the slit starts to advance requires in addition a fracture criterion (e.g. a critical
J_Ic, or a critical strain/stress over a microstructural length scale) — none is supplied —
so the load at onset of advance cannot be inferred from these results alone.
