# Answer

## 1. Requested numbers / decisions

| Quantity | Value |
|---|---|
| T | **4816.47 MPa** |
| C_rr at θ=60° | **3.9421e-3 mm^1/2** |
| C_θθ at θ=60° | **6.3318e-4 mm^1/2** |
| C_γrθ = lim √r (2ε_rθ) at θ=60° | **5.7311e-3 mm^1/2** |
| C_y(90°) = lim u_y(r,90°)/√r | **9.1390e-3 mm^1/2** |
| Screening interval on θ=0 | **[0.01000 mm, 0.06133 mm]** (non-empty) |

## 2. Model chosen and why

**Constitutive reduction (tangent/plastic branch).** As r→0 the stress grows without bound, so
σ_e ≫ σ_y and the yield cutoff factor → 1. In that limit the total-strain law becomes an
*isotropic linear* law. Writing
ε_ij = A σ_ij + B σ_kk δ_ij and comparing with the given law,

- A = (1+ν)/E + (3/2)(1/T − 1/E),   B = −ν/E − (1/2)(1/T − 1/E),

so the uniaxial compliance is A+B = 1/T, i.e. the effective tangent Young modulus is **E\* = T**,
and the effective Poisson ratio is

  **ν\* = −B E\* = 1/2 + (ν − 1/2) T/E = 0.486239** for the fitted T.

(This reproduces the bilinear uniaxial law ε = σ_y/E + (σ−σ_y)/T, whose tangent is T; the
constant σ_y/E term only produces non-singular fields and drops out of the leading 1/√r term.)

**Field.** The 2-D constraint is σ_zz = 0 (plane stress), so the leading near-tip field is the
classical mode-I plane-stress LEFM field but with the *effective* moduli (E\*=T, ν\*):

  σ_ij = (K/√(2πr)) S_ij(θ),   ε_ij = (1/T)[(1+ν\*)σ_ij − ν\* σ_kk δ_ij],

with the standard angular functions
S_xx = cos(θ/2)[1 − sin(θ/2)sin(3θ/2)], S_yy = cos(θ/2)[1 + sin(θ/2)sin(3θ/2)],
S_xy = sin(θ/2)cos(θ/2)cos(3θ/2).

**Amplitude.** Because the material is a deformation (nonlinear-elastic) theory under proportional
loading, J is path independent. For a plane-stress linear field J = K²/E\* = K²/T, so

  **K = √(J T)**.

**Treatment of ε_zz.** With σ_zz = 0, σ_kk = σ_xx+σ_yy and s_zz = −σ_kk/3, the asymptotic law gives

  ε_zz = −(ν\*/T)(σ_xx+σ_yy).

On θ=0, σ_xx=σ_yy=K/√(2πr), hence

  C_zz(0) = lim √r ε_zz = −2ν\* √(J/(2πT)).

Setting |C_zz(0)| = 0.005000 gives ν\*/√T = 0.005·√(2π)/(2√J), which solves to
**T = 4816.47 MPa** (ν\* = 0.486239, K = 62.07397 N/mm^{3/2}).

**θ = 60° coefficients.** With pref = K/√(2π) and S_ij(60°):

  C_rr = (pref/T)(S_rr − ν\* S_θθ),  C_θθ = (pref/T)(S_θθ − ν\* S_rr),
  C_γrθ = (pref/T)·2(1+ν\*) S_rθ.

**C_y(90°).** The mode-I plane-stress displacement field with (T, ν\*) gives
κ\* = (3−ν\*)/(1+ν\*), and at θ=90°:
C_y(90°) = K(3−ν\*)/(2 T √π).

**Screening interval.** On θ=0 the leading principal strain coefficients are
C_xx = C_yy = (1−ν\*)K/(T√(2π)) (positive) and C_zz = −0.005 (negative, and the largest in
magnitude). The largest absolute principal strain at radius r is therefore 0.005/√r, so the
≤0.0500 screen gives r ≥ (0.005/0.05)² = 0.01000 mm. The σ_e ≥ 10σ_y screen, σ_e = K/√(2πr),
gives r ≤ K²/(2π·(10σ_y)²) = K²/(2π·100²) = 0.06133 mm. The cap r ≤ 0.200 mm is inactive. The
interval is non-empty: **[0.01000, 0.06133] mm**.

## 3. Assumptions changed and why

- I first mis-identified the effective Young modulus as TE/(2E−T); re-deriving A+B gave E\*=T.
  The corrected value is used (confirmed by the J-integral check below).
- The 3-D law is used exactly as given; no incompressibility is imposed. ν\* < 1/2 because the
  plastic term is deviatoric only while ε_xx/ε_zz still include elastic contributions.
- For the θ=0 screening, "leading field" is read literally: the leading (asymptotic) strain
  field is used. Using the full (1−σ_y/σ_e) factor with the leading stresses instead would give
  r ≥ 0.00927 mm; the literal leading-field reading gives the clean r ≥ (0.005/0.05)² = 0.01000,
  which is what is reported.
- No T-stress / non-singular term is added; it does not affect the 1/√r coefficients.

## 4. Verification

1. **Path-independent J-integral (not a substitution into C_zz).** Using the leading stress and
   displacement fields (effective moduli T, ν\*), I evaluated
   J = ∫_Γ (W N_x − t_i u_{i,x}) ds by direct numerical contour integration on circular arcs of
   radii r = 0.01, 0.05 and 0.15 mm. All give **J = 0.800000** (to 1e-11), confirming
   E\* = T (i.e. J = K²/T) and the amplitude K = √(JT).
2. The displacement-based strain field agrees with the stress-law strains to ~1e-10 at
   θ = 10°, 60°, 90°, 150°, confirming the (T, ν\*) mode-I field is self-consistent.
3. C_y(90°) from the closed form equals a finite-difference derivative u_y/√r obtained from the
   displacement field (9.13902e-3 both ways).
4. The θ=60° polar coefficients were independently reproduced by transforming the Cartesian
   strain coefficients C_xx, C_yy, C_xy.
5. Substituting T back into 2ν\*√(J/(2πT)) reproduces 0.005000 exactly (by construction of the
   root solve, but this confirms numerical consistency).

## 5. What the screens establish; crack initiation

The three screens define an annulus (0.0100–0.0613 mm on θ=0) around the tip in which
(i) σ_e ≥ 10σ_y (the yield term is negligible, so the asymptotic tangent-field description is
appropriate), (ii) all principal strains ≤ 5% (infinitesimal kinematics are defensible), and
(iii) r ≤ 0.2 mm (small-scale-yielding match to the outer body is accurate). They are validity
checks on the leading field, not failure criteria and not error bounds for the finite-r solution.

The supplied data do **not** determine crack initiation. The constitutive law, J and C_zz give the
near-tip stress/strain state, but initiation also requires a fracture/toughness criterion — a
critical J (or K), a critical strain/distance, or a damage parameter — none of which is specified
(crack growth, damage and inertia are explicitly excluded). Hence the model can locate and
quantify the K-controlled region but cannot by itself predict when the crack advances.

## 6. Missing / uncertain

- The main ambiguity is item 4's interpretation of "leading field" (asymptotic strain vs. full
  constitutive law at the leading stress). The literal leading-field reading is reported.
- The far-field geometry (finite body, T-stress) is not given; corrections to the leading field
  are stated to be negligible for r ≤ 0.2 mm, so they do not affect the requested coefficients.
- The precision of the supplied C_zz (4 significant figures) limits T to about the same relative
  precision, well inside the 1% tolerance.
