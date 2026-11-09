# Slip and strain contours near a stationary slit — results

## Seven requested numbers

| # | Quantity | Value |
|---|----------|-------|
| 1 | Δw(r = 0.0200 mm) | **1.3181 × 10⁻³ mm** (0.0013181) |
| 2 | Δw(r = 0.0800 mm) | **1.7937 × 10⁻³ mm** (0.0017937) |
| 3 | τ_rz(r = 0.0200 mm, θ = 60°) | **+38.59 MPa** |
| 4 | τ_θz(r = 0.0200 mm, θ = 60°) | **+138.01 MPa** |
| 5 | r where \|γ\| = 0.0100, θ = 0° | **0.04689 mm** |
| 6 | r where \|γ\| = 0.0100, θ = 90° | **0.02507 mm** |
| 7 | r where \|γ\| = 0.0100, θ = 180° (face limit) | **0.01340 mm** |

All |γ| = 0.01 radii lie below 0.30 mm, inside the region where the leading field applies. Both
stresses at 60° are positive with the sign convention e_r = (cosθ, sinθ), e_θ = (−sinθ, cosθ);
w > 0 on the upper face, consistent with the stated loading.

## Governing model and how the loading fixes its magnitude

With only u_z = w(x,y) nonzero, equilibrium and the proportional power law give

- τ = τ0 (γ/γ0)^{1/n} parallel to γ (from γ = γ0(τ/τ0)^n τ/τ0),
- div τ = 0, γ = ∇w, W = ∫τ·dγ = (n/(n+1)) τ0 γ0^{−1/n} γ^{(n+1)/n} (hyperelastic, so J is path independent),
- hence div( |∇w|^{(1−n)/n} ∇w ) = 0, a p-Laplace-type equation with p = (n+1)/n = 9/7.

The leading near-tip field is the separable (HRR-type mode III) eigenfunction

  w = A r^s f(θ), γ = A r^{s−1} G(θ), G = (s²f² + f′²)^{1/2}, s = 1/(n+1) = 2/9,

where f solves the nonlinear angular ODE (obtained by substitution into the equation above)

  f″ = −f [ (1/n − 1) s² f′² + c0 (s²f² + f′²) ] / ( s²f² + f′²/n ), c0 = s(1 + (s−1)/n),

with f(0) = 0, f′(0) = 1 (antisymmetry w(x,−y) = −w(x,y); normalization absorbed into A) and
f′(π) = 0 (traction-free faces, since t_z = ∓σ_θz(±π)). This is a nonlinear eigenvalue problem
for s; shooting (DOP853 + Brent, tolerances 1e-12/1e-14) gives the single root
s* = 0.22222222222220 in (0.03, 0.6), matching the known exponent 1/(n+1) = 0.22222... to 2×10⁻¹⁴.
Numerically f(π) = 1.698433, G(0) = 1, G(90°) = 0.614353, G(180°) = 0.377430.

Magnitude from the load: for this field J = τ0 γ0^{−1/n} A^{(n+1)/n} I_ang, with

  I_ang = ∫_{−π}^{π} [ (n/(n+1)) G^{(n+1)/n} cosθ − s f G^{1/n−1} (s f cosθ − f′ sinθ) ] dθ = 2.019595,

so J = 0.120 N/mm fixes A = [ J / (τ0 γ0^{−1/n} I_ang) ]^{n/(n+1)} = 9.25609×10⁻⁴ mm^{1−s}.
This is the standard result that J alone sets the amplitude of the crack-tip field in a
deformation-theory power-law solid under a fixed loading pattern; all requested quantities follow:
Δw = 2A r^s f(π); τ_rz = τ0(A r^{s−1}G/γ0)^{1/n} s f/G, τ_θz = τ0(A r^{s−1}G/γ0)^{1/n} f′/G;
r(θ) from A r^{s−1} G(θ) = 0.01 ⇒ r = (A G/0.01)^{1/(1−s)}. Code: solve.py.

## Assumptions changed

None of the problem's assumptions was changed. Interpretations I made explicit: (i) monotonic
proportional loading ⇒ deformation theory is exact here, so W is a stored-energy density and J is
path independent; (ii) the "leading field" is the exact separable eigenfunction above, used with
the given J as amplitude — valid for r ≪ 0.30 mm; (iii) no rigid z-translation is added (w → 0 as
r → 0 already). The law has no elastic term, so no extra fitting parameter exists.

## Verification

1. **Eigenvalue vs theory:** shooting yields s* = 1/(n+1) to 2×10⁻¹⁴ (independent confirmation of
   the singularity exponent and of the angular ODE).
2. **Energy flux (path independence):** integrating W N_x − t_z w_,x over circles r = 0.02 and
   0.15 mm and squares of half-side 0.10 and 0.25 mm (contours crossing the slit line contribute
   nothing since t_z = 0 there) gives 0.120000000 N/mm in every case — the supplied J.
3. **Linear limit (n → 1):** the same pipeline returns s = 0.5, I_ang = π, f(π) = 2 exactly, and
   with μ = τ0/γ0 = 40 GPa and J = K²/(2μ) reproduces the closed-form mode-III results
   Δw = K√(8r/π)/μ (agreement to 10 digits) and K = 2J/... , i.e. the classic K-field.
4. **PDE residuals by finite differences** at interior points: |div τ| / (τ/r) ≲ 3×10⁻¹⁰ and
   |curl γ| ≲ 3×10⁻¹⁰ (roundoff level).
5. **Boundary residual:** |σ_θz|/|σ_rz| on the faces ≈ 7.5×10⁻¹⁶; |γ| evaluated at the three
   reported radii returns 0.0100000 (≤ 2.4×10⁻⁵ relative, pure print rounding).
6. Self-consistency: Δw(0.08)/Δw(0.02) = 1.360790 = 4^s exactly, as required by self-similarity.

## Missing / uncertain

- Nothing else is needed for the seven numbers (all within 2% tolerance by construction of the
  checks above); remaining uncertainty is only numerical, ≲10⁻⁹ relative, far below 2%.
- The outer geometry enters only through J = 0.120; if "negligible corrections" outside r = 0.30 mm
  were not negligible, A (hence all numbers, as A^{1/n}-type scalings) would shift by the same
  relative amount as J.
- A pure power law has no elastic range; real data below γ ≈ 0.002 would deviate, but that is
  outside the specified constitutive law.

## Do these results determine the load for crack advance?

No — the fields are known only as functions of the given energy flux J = 0.120 N/mm, and
initiation requires an additional fracture criterion (e.g., a critical J, a critical
stress/strain at a finite distance, or a cohesive-law characteristic length), none of which is
supplied.
