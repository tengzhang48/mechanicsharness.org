# Answer

## Requested numbers

| # | Quantity | Value |
|---|---|---|
| 1 | T | **4816 MPa** (4816.47 MPa) |
| 2 | C_rr (θ=60°) | **3.942×10⁻³ mm^{1/2}** |
| 2 | C_θθ (θ=60°) | **6.332×10⁻⁴ mm^{1/2}** |
| 2 | C_γrθ (θ=60°) = lim √r (2ε_rθ) | **5.731×10⁻³ mm^{1/2}** |
| 3 | C_y(90°) = lim u_y(r,90°)/√r | **9.139×10⁻³ mm^{1/2}** |
| 4 | radius interval on θ=0 passing all three screens | **[0.0100 mm, 0.0613 mm]** (non-empty) |

Derived constants used: ν_eff = 0.486239, E_eff = T, K_tip = √(J T) = 62.0740 MPa·mm^{1/2},
P = 2ν/E + 1/T − 1/E = 2.01907×10⁻⁴ MPa⁻¹.

---

## Model and constitutive reduction

**The material law is a bilinear (linear-hardening) deformation theory.** For uniaxial tension
σ: the last term of the given law gives ε = σ/E for σ ≤ σ_y and ε = σ_y/E + (σ − σ_y)/T for
σ > σ_y (verified by direct substitution: (3/2)s_xx = σ, so the plastic increment is
(1/T − 1/E)(σ − σ_y)). Thus T is the tangent (hardening) modulus and σ_y the yield stress.

**Near-tip reduction.** Loading is monotonic and proportional, mode I by symmetry, plane stress
(σ_zz = 0, ε_zz free). As r → 0 with σ = O(r^{−1/2}), σ_e → ∞, so the Macaulay factor
⟨1 − σ_y/σ_e⟩ → 1− (the −σ_y part is a bounded, non-singular strain). The leading deviatoric
response is then linear:

  ε_dev = D s,   D = (1+ν)/E + (3/2)(1/T − 1/E),
  ε_kk  = (1−2ν)/E · σ_kk   (unchanged, spherical part).

This is exactly an isotropic linear-elastic solid with

  E_eff = 3/(2D + B) = T,   ν_eff = (D − B)/(2D + B) = 0.486239,   B = (1−2ν)/E.

The step 2D + B = 3/T is exact, so **the fully-plastic near-tip tangent modulus is E_eff = T**,
not E. The associated crack field is the mode I Williams field with p = 1/2; its stress angular
functions are independent of the elastic moduli. **The stress exponent is 1/2** (linear-hardening
plastic branch ⇒ no HRR power-law exponent; the singular terms of elastic and plastic branches
scale alike).

**J, amplitude, and ε_zz.** Deformation theory has a strain-energy potential, so J is
path-independent; the near-tip J computed with the effective solid equals the supplied
J = 0.800 N/mm. In plane stress J = K_tip²/E_eff = K_tip²/T, hence

  K_tip = √(J T).                                             (a)

On θ=0, by mode-I symmetry σ_xy = 0 and the leading field has σ_xx = σ_yy = K_tip/√(2πr),
σ_zz = 0, σ_e = σ_xx. Applying the complete law (with σ_zz = 0):

  ε_zz = −(2ν/E)σ_xx − (1/T − 1/E)(σ_xx − σ_y) → −P σ_xx,
  P = 2ν/E + 1/T − 1/E = 2.01907×10⁻⁴ MPa⁻¹.                (b)

(Note 2ν_eff/T = P, so the reduced and exact expressions carry the same singular coefficient.)
Combining (a),(b):

  C_zz = lim √r ε_zz = −P K_tip/√(2π)  ⇒  P² J T = 2π C_zz².

Solving (physical root 0 < T < E; the other root is 6.36×10⁶ MPa > E):

  **(2ν/E − 1/E + 1/T)² · J · T = 2π (0.005)²  ⇒  T = 4816.47 MPa.**

## 2 and 3: strain and displacement coefficients

Mode I Williams angular stress functions S_ij(θ) = σ_ij·√(2πr)/K_tip:
S_rr = cos(θ/2)[1+sin²(θ/2)], S_θθ = cos(θ/2)[1−sin²(θ/2)], S_rθ = sin(θ/2)cos²(θ/2).
Strains from the reduced law (plane stress): ε_rr=(σ_rr−ν_eff σ_θθ)/E_eff, etc. Since
σ ~ r^{−1/2}, C = lim √r ε. At θ=60° (S_rr=1.08253, S_θθ=0.64952, S_rθ=0.375):

  C_rr  = K_tip/(E_eff√(2π))·(S_rr − ν_eff S_θθ) = 3.942×10⁻³,
  C_θθ  = K_tip/(E_eff√(2π))·(S_θθ − ν_eff S_rr) = 6.332×10⁻⁴,
  C_γrθ = 2(1+ν_eff)K_tip/(E_eff√(2π))·S_rθ       = 5.731×10⁻³  mm^{1/2}.

Plane-stress mode-I displacements, u_i = (K_tip/2G_eff)√(r/2π) f_i(θ,κ_eff),
κ_eff=(3−ν_eff)/(1+ν_eff), 2G_eff=E_eff/(1+ν_eff). At θ=90° this gives
u_y/√r = K_tip κ_eff (1/√2)/(2G_eff√(2π)) = **9.139×10⁻³ mm^{1/2}**.

## 4: radius interval on θ=0 (leading field)

- σ_e ≥ 10σ_y = 100 MPa: σ_e(θ=0)=K_tip/√(2πr) ≥ 100 ⇒ r ≤ r_a = K_tip²/(2π·100²) = **0.061325 mm**.
- largest |principal strain| ≤ 0.0500. On θ=0 the strain is diagonal; leading values are
  ε_xx=ε_yy=(1−ν_eff)σ/T and ε_zz=−2ν_eff σ/T. Since ν_eff=0.486 > 1/3, the largest is
  |ε_zz| = 2ν_eff σ/T = 0.9725 σ/T. Requiring ≤0.05 gives σ ≤ 247.64 MPa ⇒
  r ≥ r_b = K_tip²/(2π·247.64²) = **0.0100 mm**.
- r ≤ 0.200 mm.

Both r_b ≤ 0.200 and r_b ≤ r_a, so the interval is **r ∈ [0.0100, 0.0613] mm** (non-empty).
(If instead the *exact* law, without the r→0 reduction, is applied at finite r, the lower
endpoint moves to 0.00927 mm, a ~7% change; see "uncertainties".)

**What the screens establish.** They delimit a self-consistency window for the leading field: at
these radii the field is strongly yielded (σ_e ≥ 10σ_y, so the ⟨1−σ_y/σ_e⟩ factor and hence the
"fully-plastic" reduction are accurate), the strains are still small (max |ε_p| ≤ 0.05, so
infinitesimal kinematics used by the field is not violated), and r ≤ 0.200 mm where matching
corrections are stated negligible. They are checks of the field's own domain, **not** error
estimates and **not** a failure criterion.

---

## Verification

1. **Independent check suite** (`knowledge/checks/crack_tip_checks.py`, `verify.py`): the mode-I
   field with p=1/2, law="linear", nu=ν_eff, plane="stress" PASSES every applicable check —
   traction-free faces, mode-I symmetry, angular and Cartesian equilibrium, Saint-Venant
   compatibility, finite energy, crack exponent fixed by J, crack faces open, and the
   displacement/shear-relation compatibility. This confirms the field solves the stated model
   with the reduced law.
2. **J numerically, not by substitution** (`solve.py`): integrating J = ∫(W N_x − t·u_x)ds on
   circles r = 0.1, 1, 10 with the effective solid returns 0.8000000000 every time, equal to
   K_tip²/T. This validates the J ⇔ K ⇔ T relation and the displacement field (including κ_eff).
3. **Strain coefficients by finite differences** (`verify2.py`): numerical derivatives of the
   displacement field at θ=60° reproduce C_rr, C_θθ, C_γrθ to 9 digits; C_y(90°) likewise.
4. **Reduction checked against the full law**: at σ=10⁷ MPa (where the yield correction is
   negligible) the full law gives ε_zz/ε_xx = −1.89285794, matching the reduced value
   −2ν_eff/(1−ν_eff) = −1.89285804.
5. The given datum is reproduced exactly: the solved T returns C_zz = −0.0050000000.

## What is missing or uncertain

- **Part 4 lower endpoint.** "Evaluated using the leading field" is read as the asymptotic
  (reduced) field, giving r_b = 0.0100 mm. If one instead inserts the leading *stress* into the
  *unreduced* complete law at finite r, r_b = 0.00927 mm (the two differ because the r→0 limit
  drops the bounded −σ_y term). The coefficients in parts 2–3 are unambiguous (they are defined
  by the r→0 limit, which forces the reduced law).
- **Applicability.** The leading field is an intermediate asymptotic: it presumes small-scale
  yielding, proportional monotonic loading (so deformation theory ≈ flow theory), J-dominance,
  and radii large compared with the fracture-process/blunting zone. Here δ_t ~ d_n J/σ_y with
  d_n ≈ 0.3–0.8, i.e. δ_t ≈ 0.02–0.06 mm, which is *larger* than r_b = 0.01 mm; the inner screen
  endpoint may therefore already sit inside the blunting/process zone where small-strain
  kinematics is questionable. None of this is established by the supplied data; the problem
  asserts the leading field is adequate for r ≤ 0.200 mm.
- **The data do not determine crack initiation.** T, J, and C_zz are field/amplitude parameters.
  A fracture criterion needs a toughness (critical J or K) or a coupled stress+energy condition
  for the notch/process zone; no such data are given. The screens of part 4 are consistency
  checks on the asymptotic field, not initiation criteria.

## Files
- `solve.py` — main computation.
- `verify.py` — runs the provided constraint checks.
- `verify2.py` — finite-difference and full-law cross-checks.
