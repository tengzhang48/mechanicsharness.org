# Infer T from a transverse strain field — answer

Units: mm, N, MPa. Strains are tensor strains; engineering shear = 2 x tensor shear.

## Requested numbers

| # | Quantity | Value |
|---|---|---|
| 1 | T | **4816 MPa** (4816.472) |
| 2 | C_rr(60°) = lim √r ε_rr | **3.942e-3 mm^1/2** |
| 2 | C_θθ(60°) = lim √r ε_θθ | **6.332e-4 mm^1/2** |
| 2 | C_γrθ(60°) = lim √r (2 ε_rθ) | **5.731e-3 mm^1/2** |
| 3 | C_y(90°) = lim u_y(r,90°)/√r | **9.139e-3 mm^1/2** |
| 4 | radii on θ=0 passing all three screens | **[0.01000, 0.06133] mm** (non-empty) |

Intermediate: K = √(J T) = 62.074 MPa·mm^1/2.

## Model and why

**Asymptotic reduction of the total-strain law.** Near the tip the stresses are singular, so
σ_e → ∞ and ⟨1 − σ_y/σ_e⟩ → 1. The last term therefore becomes linear in s_ij and the
*complete* law reduces, at leading order, to an anisotropic-in-shear linear elastic law:

  ε_ij = (1+ν)/E σ_ij − (ν/E) σ_kk δ_ij + (3/2)(1/T − 1/E) s_ij.

With plane stress (σ_zz = σ_xz = σ_yz = 0; ε_zz free) this is standard 2D isotropic plane
stress in disguise. Let S = σ_xx + σ_yy and

  α1 = (1+ν)/E + (3/2)(1/T − 1/E),   α2 = ν/E + (1/2)(1/T − 1/E).

Then, for in-plane a,b and using s_ab = σ_ab − S/3 δ_ab,

  ε_ab = α1 σ_ab − α2 S δ_ab,   ε_zz = −α2 S.

This is the 2D law of a standard isotropic plane-stress material with effective in-plane
moduli (verified below)

  E2 = 1/(α1 − α2) = T,     ν2 = α2 T = 1/2 + (T/E)(ν − 1/2),
  ε_zz = −ν2 S / E2 = −α2 S   (matches the ε_zz above exactly).

**Field.** Because the asymptotic material is homogeneous, isotropic and linear, the leading
mode-I field is the classical Williams field with p = 1/2 and angular functions independent of
the moduli (they follow from equilibrium, compatibility and the traction-free faces only):

  σ_ij = (K/√(2π r)) f_ij(θ),
  f_xx = cos(θ/2)[1 − sin(θ/2)sin(3θ/2)],
  f_yy = cos(θ/2)[1 + sin(θ/2)sin(3θ/2)],
  f_xy = sin(θ/2)cos(θ/2)cos(3θ/2).

The plane-strain HRR machinery of the supplied knowledge is *not* used: the extra term is
linear in stress, so it does not dominate the elastic term or change the exponent, it only
changes the effective compliance. (J2 pure-power hardening would instead give p = 1/(n+1).)

## Fixing T

On θ = 0 the field is equibiaxial in-plane with σ_xx = σ_yy = K/√(2πr), σ_zz = 0, so
S = 2K/√(2πr) and ε_zz = −2α2 K/√(2πr). Then

  C_zz(0) = −2 α2 K/√(2π) = −0.0050.

For a straight crack in this homogeneous linear material J = K²/E2 = K²/T (verified by direct
contour integration below), and J is path independent for the deformation-theory law, so the
supplied J applies at the tip:

  J = K²/T = 0.800.

Solving the two equations (α2 = 1/(2T) + (ν−1/2)/E):

  T = 4816.472 MPa,   K = √(0.8·T) = 62.074 MPa·mm^1/2.

## The other coefficients

At θ = 60°: f_rr = 1.082532, f_θθ = 0.649519, f_rθ = 0.375000, S_ang = f_rr+f_θθ = 1.732051,
and α1 = 3.085741e-4, α2 = 1.009533e-4, K/√(2π) = 24.7640:

  C_rr = (K/√(2π))(α1 f_rr − α2 S_ang)  = 3.942e-3,
  C_θθ = (K/√(2π))(α1 f_θθ − α2 S_ang)  = 6.332e-4,
  C_γrθ = (K/√(2π))(2 α1 f_rθ)          = 5.731e-3.

For C_y(90°): on θ = 90° e_r = ŷ, so u_y = u_r, and the leading field has u_r = 2 C_rr r^1/2
(since ε_rr = ∂u_r/∂r gives C_rr = U_r/2). Numerically C_rr(90°) = 4.5695e-3 and

  C_y(90°) = 2 C_rr(90°) = 9.139e-3.

This equals the standard plane-stress solution u_y = (K/(2μ2))√(r/2π)·(√2/2)κ2 with
μ2 = T/(2(1+ν2)) = 1620.36, κ2 = (3−ν2)/(1+ν2) = 1.69136, ν2 = 0.486239.

## Screens on θ = 0

Leading field: σ_xx = σ_yy = σ(r) = K/√(2πr), σ_zz = 0.

- σ_e at θ = 0 equals σ(r) (equibiaxial in-plane, zero through-thickness stress), so
  σ_e ≥ 10σ_y = 100 MPa gives r ≤ (K/√(2π)/100)² = **0.061325 mm**.
- Principal strains are {ε_xx, ε_yy, ε_zz} = {(α1−2α2)σ, (α1−2α2)σ, −2α2 σ}; the largest
  magnitude is |ε_zz| = 2α2 σ. The measured C_zz fixes 2α2 K/√(2π) = 0.005, so
  |ε_zz| = 0.005/√r. The condition max|ε_principal| ≤ 0.0500 gives √r ≥ 0.1, i.e.
  r ≥ **0.010000 mm**.
- r ≤ 0.200 mm is looser than the σ_e bound.

Intersection: **[0.01000, 0.06133] mm**, non-empty.

## Verification (not just substitution into C_zz)

1. **Direct J contour integral.** Reconstructing the near-tip field (stresses and the
   effective plane-stress strains) and integrating J = ∫_Γ (W N_x − t_i u_{i,x}) ds on
   circles R = 0.1, 0.5, 1.0, 2.0 mm (Cartesian ∂u/∂x, so no dropped polar-basis terms) gives
   **0.800000** on every contour, matching the supplied J. This checks the independent
   J = K²/T relation, not the C_zz data.
2. **Topic checks** (`knowledge/checks/crack_tip_checks.py`, `check_field` with
   law="linear", plane="stress", ν = ν2, p = 0.5): traction-free faces (3e-23), mode-I
   symmetry (1e-16), angular equilibrium (1e-10), Cartesian equilibrium (1e-8),
   Saint-Venant compatibility (4e-7), finite energy, crack exponent p = 1/2, crack faces
   open, and displacement/shear compatibility (3e-12) all PASS.
3. **Cross-check of C_y**: the analytic mode-I displacement formula and the strain-integration
   route (u_y = u_r = 2 C_rr r^1/2) agree to 1e-15.
4. Identity α1 − α2 = 1/T and ν2 = α2 T reproduce E2 = T and ε_zz = −ν2 S/E2, confirming the
   effective-modulus identification.

## Assumption changed and why

- *Changed:* the full bilinear total-strain law is replaced by its σ_e → ∞ linear limit at the
  tip. This is a controlled asymptotic step (the dropped factor is 1 − σ_y/σ_e → 0 relative to
  the retained term), and it is what makes p = 1/2 and a K-field correct here.
- *Not changed:* plane stress (σ_zz = 0, ε_zz free) and infinitesimal kinematics, as specified.

## What the screens establish; crack initiation

The screens are consistency windows, not error bars: they locate radii where (i) the field is
in a J-dominated, high-triaxiality plastic regime far above yield, (ii) strains remain within
the small-strain assumption, and (iii) the field is inside the range where matching to the
surrounding body keeps leading-term corrections negligible. They say nothing about failure.

The supplied constitutive law and loading (J) fix the amplitude K and the entire near-tip
field, but **do not determine crack initiation**: no fracture criterion (e.g. a critical J,
CTOD, or a process-zone/cohesive strength) and no material toughness are given, and J is only
a loading/field parameter. Deciding initiation requires additional data beyond E, ν, σ_y, T
and J.

## Uncertainties / missing

- The result relies on taking the σ_e → ∞ limit; at the σ_e = 10σ_y screen boundary the
  neglected correction to the law is ~10%, so the leading field is used as the problem
  defines it.
- The identification J = K²/T assumes the standard plane-stress energy-release relation for
  the effective in-plane material; it was verified numerically (contour integral) but is an
  asymptotic statement.
- Deformation theory vs. flow theory and any crack-tip blunting/process zone are outside the
  model; no initiation threshold is inferred.
