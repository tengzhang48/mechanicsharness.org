# Slip and strain contours near a stationary slit — answer

## 1. Requested numbers (mm, N, MPa; angles as stated)

| # | Quantity | Value |
|---|---|---|
| 1 | Δw(r = 0.0200 mm) = w(r,+π) − w(r,−π) | **1.3181×10⁻³ mm** |
| 1 | Δw(r = 0.0800 mm) | **1.7937×10⁻³ mm** |
| 2 | τ_rz(r = 0.0200 mm, θ = 60°) | **+38.59 MPa** |
| 2 | τ_θz(r = 0.0200 mm, θ = 60°) | **+138.01 MPa** |
| 3 | r where |γ| = 0.0100, θ = 0° | **4.689×10⁻² mm** |
| 3 | r where |γ| = 0.0100, θ = 90° | **2.507×10⁻² mm** |
| 3 | r where |γ| = 0.0100, upper face θ = 180° | **1.340×10⁻² mm** |

(Both stresses are positive: with the upper face driven in +z, the shear vector points
outward/upward in the upper half-plane.)

## 2. Model and how the loading fixes its magnitude

This is an anti-plane (mode III) crack problem for a pure power-law deformation-plasticity
material, so I used the separable leading near-tip field (the anti-plane analogue of the
HRR field):

  w(r,θ) = a r^p h(θ),  p = 1/(n+1) = 0.22222 (n = 3.5).

Why this model:

* Equilibrium div τ = 0, compatibility γ = ∇w, and the given law
  τ_vec = τ0(|γ|/γ0)^{1/n} γ_vec/|γ| admit a separable solution with
  τ ∝ r^{−1/(n+1)}, γ ∝ r^{−n/(n+1)}. The exponent p = 1/(n+1) is fixed by requiring a
  finite, nonzero, path-independent J (energy density W ∝ r^{−1}); this is the standard
  crack result, also stated in `knowledge/TOPIC.md` §3. (Admissibility: energy integrable.)
* With q(θ) = √(p²h² + h′²), equilibrium reduces to
  d/dθ[q^ν h′] + β q^ν h = 0, ν = (1−n)/n = −5/7, β = n/(n+1)² = 0.17284,
  with h(0) = 0 (antisymmetry) and h′(π) = 0 (traction-free face). The Eshelby
  conservation law behind J gives a pointwise first integral,
  sinθ h′² − cosθ h h′ − β sinθ h² = 0, i.e. the first-order ODE
  h′ = h[cosθ + √(cos²θ + 4β sin²θ)]/(2 sinθ), which I integrated (DOP853, tol 1e-12),
  normalized h ~ θ at θ = 0. Both boundary conditions are then automatic.
* Amplitude from the applied energy flux: evaluating J on a circular arc gives
  J = τ0γ0 (a/γ0)^{(n+1)/n} I_J with, in this normalization,
  I_J = ∫₋π^{π}[(n/(n+1)) q^{(n+1)/n} cosθ − q^ν p h(p h cosθ − h′ sinθ)]dθ = 2.0195953.
  J = 0.120 N/mm ⇒ K_H = J/I_J = 0.0594178 N/mm ⇒
  a = γ0(K_H/(τ0γ0))^{n/(n+1)} = 9.25609×10⁻⁴ mm^{7/9}.  Key field values:
  h(π) = 1.6984332; at θ = 60°: h = 0.960852, h′ = 0.763695, q = 0.792983.

The three requested items follow from Δw = 2 a r^p h(π);
τ_rz = τ0(a/γ0)^{1/n} r^{−1/(n+1)} q^ν p h (same with h′ for τ_θz);
and r = (a q(θ)/0.01)^{(n+1)/n} for the strain contour.

## 3. Assumptions changed

None. I used exactly the specified pure power law (no elastic term, as stated), the leading
local field as instructed (every requested position and contour radius — largest is
r = 0.0469 mm — lies well inside the stated 0.30 mm region where outer-geometry corrections
are negligible), and quasistatic monotonic loading, which makes deformation theory exact
within the separable zone (the field is exactly proportional at fixed θ).

## 4. Verification (all run by `solve.py` and `check_shoot.py`)

* **Known limiting solution (n = 1):** the same code with n = 1 reproduces the linear-elastic
  K_III field: h = 2 sin(θ/2) to 1.5×10⁻¹², I_J = π to 1.2×10⁻¹³, and the amplitude satisfies
  J = K_III²/(2μ), μ = τ0/γ0, to 6×10⁻¹⁴.
* **Energy-flux calculation (independent):** rebuilding w(x,y) on circles
  r = 0.02/0.05/0.10 mm, forming w_,x, w_,y by Cartesian finite differences, τ from the
  vector law and W = (n/(n+1)) τ·γ, gives J = 0.120000 N/mm (relative error ≈ 1×10⁻¹¹),
  confirming both the amplitude inversion and path independence of the stated J integral.
* **Equilibrium:** residual of the second-order ODE ≤ 4×10⁻⁹; independently shooting that
  ODE from h(0)=0, h′(0)=1 lands on h′(π) = −3×10⁻¹⁴ (i.e. p = 1/(n+1) is the correct
  eigenvalue) and matches the first-integral solution to 3×10⁻¹².
* **Boundary residual:** τ_θz = 0 on both faces (h′(π) = 0 to machine precision), and the
  antisymmetry h(−θ) = −h(θ) is built in.

## 5. Missing or uncertain

* The results are the leading singular term only; higher-order terms from the outer geometry
  are stated to be negligible for r ≤ 0.30 mm but are unquantified here.
* Small-strain theory fails within a few crack-opening displacements of the tip
  (δ_t ~ J/τ0 ≈ 1.5×10⁻³ mm); the smallest requested radius (0.02 mm) is ~13 δ_t, so the
  field is self-consistent there, but anything much closer to the tip would not be.
* **Crack advance:** No — these results alone do not determine the load at which the slit
  starts to advance, because the analysis only maps the applied energy flux J onto near-tip
  fields, while onset of advance requires an independent material fracture criterion (e.g. a
  critical J_c or a critical crack-tip slip) that is not specified.
