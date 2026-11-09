# ANSWER — infer T from a transverse strain field

## 1. Requested numbers

| Quantity | Value |
|---|---|
| **T** | **4816.5 MPa** (4800 would round to; full: 4816.47 MPa) |
| C_rr(60°) = lim √r ε_rr | **0.003942 mm^1/2** |
| C_θθ(60°) = lim √r ε_θθ | **0.000633 mm^1/2** |
| C_γrθ(60°) = lim √r (2ε_rθ) | **0.005731 mm^1/2** |
| C_y(90°) = lim u_y(r,90°)/√r | **0.009139 mm^1/2** |
| Screen interval on θ = 0 | **[0.0100 mm, 0.0613 mm]**, nonempty |

Auxiliary: near-tip effective constants E\* = T = 4816.47 MPa, ν\* = 0.486239; local stress
intensity K_loc = √(J·T) = 62.074 MPa·mm^1/2.

## 2. Model chosen and why

**Constitutive reduction.** The given law is J2 deformation (total-strain) theory with bilinear
uniaxial response: σ/E up to σ_y, then tangent T (checked: at σ = 2σ_y it returns
σ_y/E + σ_y/T exactly). Since stresses are singular like r^(−1/2), σ_e → ∞ as r → 0, so the
Macaulay factor ⟨1 − σ_y/σ_e⟩ → 1 and the omitted −σ_y/σ_e part contributes only bounded
(O(r⁰)) strain corrections. The leading field therefore obeys an isotropic linear law

  ε_ij = (1+ν\*)/E\* σ_ij − (ν\*/E\*) σ_kk δ_ij,  with **E\* = T** exactly and
  ν\* = 1/2 − (T/2E)(1 − 2ν) = 0.486239.

E\* = T is exact because the plastic strain is deviatoric (incompressible): ε_kk = (1−2ν)/E σ_kk
is purely elastic, and comparing 1/E\* = (1+ν\*)/E\* − ν\*/E\* gives 1/E\* = 1/T.

**Field.** Plane stress (σ_zz = 0, ε_zz free), stationary crack, mode I, faces traction-free,
deformation theory (strain energy exists, loading stated monotonic + proportional), so J is
path-independent and a finite nonzero J forces the crack exponent p = 1/2 (knowledge TOPIC §3).
The leading field is thus the Williams mode-I field with (E\*, ν\*): σ_ij = K_loc r^(−1/2) f_ij(θ)/√(2π),
and J = K_loc²/E\* = K_loc²/T (plane stress). This is consistent with Hutchinson 1968 eq. (42):
for a piecewise-linear material the plane-stress amplitude carries the factor √(E_tan/E), i.e.
K_loc = √(J E)·√(T/E).

**Determining T.** On θ = 0, σ_xx = σ_yy = K_loc/√(2πr), so σ_e = σ (plane stress) and
ε_zz = −(ν\*/T)(σ_xx+σ_yy) ⇒ √r ε_zz = −(2ν\*/T)K_loc/√(2π) = −0.005. With K_loc = √(0.8T) this
is one equation in T; solved by Brent root-finding: **T = 4816.47 MPa** (0 < T < E ✓).

**Treatment of ε_zz.** ε_zz is free (plane stress), and with the ⟨⟩ → 1 limit,
ε_zz = −(ν/E + ½(1/T − 1/E))(σ_xx+σ_yy) = −(ν\*/T)(σ_xx+σ_yy); generally
C_zz(θ) = −(2ν\*/T)(K_loc/√(2π)) cos(θ/2). The elastic and plastic parts both enter; only the
volumetric strain is purely elastic.

**Part 2** from σ_rr = 1.0825, σ_θθ = 0.6495, σ_rθ = 0.375 (× K_loc/√(2πr), θ = 60°) and
plane-stress Hooke with (T, ν\*): ε_rr = (σ_rr−ν\*σ_θθ)/T, ε_θθ = (σ_θθ−ν\*σ_rr)/T, 2ε_rθ = 2(1+ν\*)σ_rθ/T.

**Part 3**: at θ = 90°, u_y = u_r; with u ~ r^1/2 U(θ), U_r = 2ε_rr^ang, so
C_y(90°) = 2 C_rr(90°) = (3−ν\*)/√2 · K_loc/(T√(2π)) = 0.009139 (no rigid translation; the
symmetry u_θ(0) = 0 fixes the integration constant). Cross-checked against the closed-form
plane-stress u_r = (2K/E)√(r/2π) cos(θ/2)[(1−ν)+(1+ν)sin²(θ/2)] — identical to 6 digits.

**Part 4 screens** (θ = 0, leading field; σ = K_loc/√(2πr), principal strains
ε_1 = ε_2 = (1−ν\*)σ/T, ε_3 = ε_zz = −2ν\*σ/T with 2ν\* > 1−ν\*):
- σ_e = σ ≥ 10σ_y = 100 MPa ⇒ r ≤ K_loc²/(2π(10σ_y)²) = 0.061325 mm;
- max |principal strain| = |ε_zz| = |C_zz(0)|/√r = 0.005/√r ≤ 0.0500 ⇒ r ≥ 0.010000 mm;
- r ≤ 0.200 mm (given matching statement).
All three hold on **[0.0100, 0.0613] mm** (at r_lo: strain = 0.0500, σ_e = 247.6 MPa; at
r_hi: strain = 0.0202, σ_e = 100.0 MPa).

## 3. Assumptions changed

None of the stated physical assumptions was changed. Two working reductions, both justified:
(i) I replaced the Macaulay factor by 1 for the *leading* r → 0 coefficient — this is not an
assumption change but the limit that defines C_zz; quantitatively, the full-law vs
effective-law strain gap at θ = 60° is 3.3%, 5.7%, 8.1% at r = 0.010, 0.030, 0.061 mm,
tracking σ_y/σ_e as it should (the ⟨⟩ correction is bounded, O(r⁰) in √r-normalized units).
(ii) For J I used W = ½σ_ijε_ij of the tangent-linear near-tip field, consistent with
deformation theory since the contour can be shrunk into the fully plastic region where the law
is linear; the σ_y/σ_e cross-terms vanish as r → 0. Also, the provided checks library sets all
material constants to 1, so for its reported-J comparison I rescaled my field by 1/√T (a
normalization convention only, not a model change).

## 4. Verification

1. **Provided checks** (`knowledge/checks/crack_tip_checks.py`) on my angular field,
   law="linear", ν = ν\*, p = 0.5, plane stress, mode I, β = π: **all 10 applicable checks
   PASS** — traction-free faces (2.6e-23), mode-I symmetry (1.4e-16), angular equilibrium
   (1e-10, also tests the exponent), Cartesian equilibrium (1.1e-8), Saint-Venant
   compatibility (4.1e-7), finite energy, exponent p = 1/2 fixed by J (exact), crack faces open,
   and displacement/strain compatibility (2.6e-12).
2. **Quantitative, not substitution into C_zz**: the field's own J-integral. The checks compute
   J from my field in Cartesian components: **0.800000 on contours r = 0.5, 1, 2, matching the
   supplied J = 0.800 to 8×10⁻¹⁴**. I also coded an **independent** physical-unit contour
   integral (own strains→displacements integration, W = ½σ:ε, t_i ∂u_i/∂x by finite
   differences): **J = 0.800000 on r = 0.005, 0.02, 0.05, 0.1 mm**. (An initial wrong value
   0.5486 traced to a sign error in my polar→Cartesian stress transform, and an earlier one to
   ∂u_y/∂y instead of ∂u_y/∂x — both found by this cross-check and fixed.)
3. Uniaxial reduction of the given law is exactly bilinear (σ_y/E + (σ−σ_y)/T).
4. C_y(90°) computed two independent ways (kinematic integration identity 2C_rr(90°) and the
   closed-form u_r) agree to 6 digits (0.00913902 vs 0.009139).
5. Internal consistency of the fit: with T = 4816.47, K_loc = √(JT) reproduces the measured
   C_zz(0) = −0.005000; angular form predicts C_zz(60°) = −0.004330 (not used in the fit,
   available for external comparison).

## 5. What the screens establish; what is missing or uncertain

**Screens.** They bracket the annulus where the *leading field* is self-consistent, each row
mapping to a TOPIC §6 condition: σ_e ≥ 10σ_y puts points deep in the plastic branch (the
⟨σ_y/σ_e⟩ ≤ 0.1 correction to the linear-tangent reduction, and uniformity of the plastic
state); |principal strain| ≤ 0.05 keeps infinitesimal kinematics defensible; r ≤ 0.200 mm keeps
outer matching corrections negligible. They are defined checks on the asymptotic model, **not**
error estimates of the full finite-r solution, and they say nothing about adequacy of
small-strain continuum at the tip itself (blunting/process zone, δ_t ~ J/σ_y scale).

**Crack initiation?** No. The supplied data (E, ν, σ_y, T, J, C_zz) fix the *field* and its
amplitude, but initiation requires a failure criterion — a critical K, J, CTOD or a damage
model, i.e., a measured toughness — which is not supplied (TOPIC §6 item 9: an amplitude is a
field parameter, not a fracture criterion). Additional open points: deformation theory vs flow
theory is legitimate only under the stated proportional monotonic loading (Hutchinson §8 notes
stress proportions differ across the elastic–plastic boundary; the problem statement stipulates
proportionality, which I relied on); plane stress is a 2D idealization whose zone must exceed
the thickness; and the screens above assume the leading field alone, so near the interval
endpoints the neglected O(r⁰) terms are of the size quantified in §3(i).

Uncertainty of the answer itself: T enters through √T (from J) and 1/T (from ν\*/T), so the
solution is well-conditioned; the reported digits are far inside the 1%/2% tolerances.
