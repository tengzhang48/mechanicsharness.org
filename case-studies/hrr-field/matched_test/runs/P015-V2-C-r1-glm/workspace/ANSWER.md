# ANSWER

Units: mm, N, MPa. Given E = 70,000, ν = 0.30, σ_y = 10.0, J = 0.800 N/mm, C_zz(0) = −0.005 mm^1/2.

## 1. Results (requested numbers)

| Quantity | Value |
|---|---|
| **T** | **4816.5 MPa** (solve gives 4816.472) |
| C_rr(60°) | **3.9420 × 10⁻³ mm^1/2** |
| C_θθ(60°) | **6.3318 × 10⁻⁴ mm^1/2** |
| C_γrθ(60°) = lim √r·(2ε_rθ) | **5.7311 × 10⁻³ mm^1/2** |
| C_y(90°) = lim u_y(r,90°)/√r | **9.1390 × 10⁻³ mm^1/2** |
| Screening interval on θ = 0 | **r ∈ [0.01000, 0.06133] mm** (nonempty) |

Useful derived numbers: tangent modulus pair E_t = T = 4816.5 MPa, ν_t = 0.486239; stress intensity K_t = √(J·E_t) = 62.074 MPa·mm^1/2; K_t/(√(2π)E_t) = 5.1415×10⁻³ mm^1/2.

## 2. Model chosen and why

**Constitutive reduction (tangent moduli).** The total-strain law is a Hencky (deformation-theory, nonlinear-elastic) law; the stated monotonic, proportional loading makes it legitimate. As r → 0, σ_e → ∞, so ⟨1 − σ_y/σ_e⟩ → 1 − σ_y/σ_e and

ε_ij = (1+ν)/E σ_ij − ν/E σ_kk δ_ij + (3/2)(1/T − 1/E) s_ij − (3/2)(1/T − 1/E)(σ_y/σ_e) s_ij.

The first three terms are **linear** in σ and isotropic: ε = A σ + B σ_kk δ with A = (ν−½)/E + 3/(2T), B = (½−ν)/E − 1/(2T). Rewriting ε = (1+ν_t)/E_t σ − ν_t/E_t σ_kk δ gives exactly

- **E_t = T** (uniaxial asymptotic compliance 1/E + (1/T − 1/E) = 1/T), and
- **ν_t = ½ − (½−ν)T/E** (plastic deviatoric slope is incompressible; ν_t = ν when T = E).

The last term is bounded (s_ij/σ_e is O(1), coefficient σ_y), i.e. an r⁰ strain "shift" that does not affect any lim √r · ε coefficient, contributes O(r) displacements, and vanishes from J as r → 0 (verified numerically below).

**Field.** Because the operator is asymptotically the isotropic linear one, the dominant Williams eigenfield is the r^(−1/2) plane-stress **mode I** K-field (pure mode I by symmetry of the loading), with moduli (E_t, ν_t): σ_xx,σ_yy,σ_xy = (K_t/√2πr)·{c(1∓s·s₃), c·s·c₃}, c = cos θ/2, s = sin θ/2, s₃ = sin 3θ/2, c₃ = cos 3θ/2. The amplitude follows from path-independence of J (exact for this deformation theory): evaluated on a vanishing contour with the tangent field, **J = K_t²/E_t** (plane stress), so K_t = √(J T).

**Treatment of ε_zz.** Plane stress (σ_zz = σ_xz = σ_yz = 0) with ε_zz free: from the law, ε_zz = −(ν_t/E_t)(σ_xx + σ_yy) + O(r⁰). On θ = 0, σ_xx = σ_yy = K_t/√2πr, σ_xy = 0, σ_e = K_t/√2πr, hence

  √r ε_zz(0) = −(ν_t/E_t)·2K_t/√(2π) = −2[1/(2T) − (½−ν)/E]·√(JT)/√(2π) = −0.005.

The left side is strictly decreasing in T on (0, E) ⇒ unique root: **T = 4816.47 MPa**. (Consistency check: the full nonlinear law evaluated at θ = 0 gives √r·ε_zz → −0.004807, −0.004981, −0.0049998 at r = 10⁻², 10⁻⁴, 10⁻⁸ — the r^(−1/2) tangent limit is the measured −0.005.)

**Q2.** Polar stress coefficients at 60° (×K_t/√2πr): σ_rr = c(1+s²) = 1.082532, σ_θθ = c³ = 0.649519, σ_rθ = sc² = 0.375. Strains with (E_t, ν_t) plane stress:
C_rr = P(σ_rr − ν_tσ_θθ), C_θθ = P(σ_θθ − ν_tσ_rr), C_γrθ = 2(1+ν_t)Pσ_rθ, with P = K_t/(√2π E_t) = 5.14151×10⁻³ mm^1/2 → the values in the table.

**Q3.** Integrating the strain field (u_y(0,y) = ∫₀^y ε_yy dy′, u_y(0,0) = 0) gives exactly C_y = 2C_yy(90°) = (3−ν_t)·P·sin 45° = 9.1390×10⁻³ mm^1/2. (I also constructed the full displacement field from the polar strain relations, u = P√r[2g_rr e_r + 4(g_rr′−g_rθ)e_θ]; its u_y(r,90°)/√r and its CTOD, 8K_t/E_t√(r/2π), match the textbook plane-stress values.)

**Q4.** On θ = 0 the leading field gives σ_e = K_t/√2πr, principal strains ε₁ = ε₂ = (1−ν_t)K_t/(E_t√2πr), ε₃ = ε_zz = −2ν_t K_t/(E_t√2πr); since 2ν_t = 0.972 > 1−ν_t = 0.514, the largest |principal strain| is |ε_zz| whose √r-coefficient is exactly |C_zz| = 0.005. So:
- σ_e ≥ 10σ_y = 100 MPa ⟺ r ≤ K_t²/(2π·10⁴) = JT/(2π(10σ_y)²) = 0.061325 mm;
- |ε|_max ≤ 0.05 ⟺ r ≥ (0.005/0.05)² = 0.010000 mm;
- r ≤ 0.200 mm non-binding.
Intersection: **[0.010000, 0.061325] mm**, nonempty (at r_hi σ_e = 100.0 MPa exactly; at r_lo |ε_zz| = 0.0500).

## 3. Assumptions I changed, and why

- I first reduced the asymptotic law to a *secant-like* effective modulus E_eff = (1/E + 1/T)⁻¹; re-expanding the law showed the elastic part of the compliance cancels in the slope and the correct asymptotic (tangent) pair is **E_t = T, ν_t = ½ − (½−ν)T/E**. This moved T from ≈5175 to ≈4816 MPa. Confirmed by feeding the leading stress field through the *complete* law and checking the r → 0 limits.
- I first used a displacement formula recalled from memory, (1+ν)(K/E)√(r/2π)·[cos(θ/2)(κ−1+2sin²(θ/2)), sin(θ/2)(κ+1−2cos²(θ/2))]. It passed strain checks only at shear-free rays and gave J = 0.0885 instead of 0.800 on the contour. I discarded it and reconstructed u from the polar strain–displacement relations (u = P√r[2g_rr e_r + 4(g_rr′−g_rθ)e_θ]), which matches the constitutive strains at all θ and reproduces CTOD = 8K_t/E_t√(r/2π) exactly. Q1/Q2/Q4 never used displacements; Q3 is anchored on the y-axis identity C_y = 2C_yy(90°), which is independent of the rotation field.

## 4. Verification (beyond substituting into C_zz)

1. **J-contour integral of the full nonlinear law** (code: `verify3.py`): integrating J = ∮(W N_x − t_i u_{i,x})ds on circles of radius ρ with W the exact work density of the *complete* law along the proportional path (closed form W = ½σ:ε^e + (3/4)(1/T−1/E) s:s·(1−(σ_y/σ_e)²)) gives 0.799835, 0.799995, 0.800000, 0.800000 at ρ = 10⁻¹, 10⁻², 10⁻³, 10⁻⁴ mm (error ∝ √ρ → 0), and exactly 0.800000 with the tangent-modulus W at every ρ. This verifies the field, its amplitude K_t = √(JT), and the plane-stress relation J = K_t²/E_t — independent of the C_zz calibration.
2. **Full-law limits at θ = 60°** (code: `solve.py`): the complete constitutive law evaluated on the leading stresses gives √r·(ε_rr, ε_θθ, 2ε_rθ) → (3.9419, 0.6332, 5.7309)×10⁻³ at r = 10⁻⁸ mm, converging to the reported coefficients — checks the polar transformation and the neglect of the bounded terms.
3. **Displacement field**: strains from the reconstructed u match the law at θ = 0, 30, 45, 60, 90, 120° to 7+ digits; u_x(r,±π) = 0 (no interpenetration), u_y > 0 (opening), CTOD = 8K_t/E_t·√(r/2π) exactly; σ_θθ = σ_rθ = 0 on the faces (traction free); C_y = φ₀F_y(90°) = 2C_yy(90°) by two independent routes.
4. **Uniqueness of T**: the calibration function is strictly decreasing on (0, E) (checked numerically), so the root is unique.

## 5. What the screens establish; is initiation determined?

The three screens are defined self-consistency checks for using the leading field, not error estimates: σ_e ≥ 10σ_y keeps σ_y/σ_e ≤ 0.1, i.e. the neglected bounded terms are ≤ ~10% of the strain (deep plastic asymptote); |ε|_max ≤ 0.05 keeps infinitesimal kinematics defensible; r ≤ 0.200 mm is where matching to the surrounding body is stated negligible. The interval [0.0100, 0.0613] mm is the region where all three hold simultaneously — e.g. below 0.0100 mm the model strains exceed 5% and small-strain theory loses meaning (the formal r → 0 limit is not physically accessible).

**Crack initiation is NOT determined** by the supplied data. The constitutive law is damage-free and contains no failure criterion (no K_c, J_c, CTOD_c, critical strain, or cohesive/damage parameters); J = 0.800 N/mm is a loading (amplitude), not a resistance. Also, deformation theory is valid only under the monotonic proportional loading assumed here, while initiation/growth involves unloading and path changes. One would need an additional initiation/toughness criterion to say anything about initiation.

## 6. Missing / uncertain

- Nothing in the requested outputs is uncertain beyond the stated tolerances; T is determined by a monotone scalar equation (4816.5 MPa; 1% band ≈ 4768–4865).
- J is trusted as supplied ("external calculation"); the negligibility of matching corrections for r ≤ 0.200 mm is taken as given.
- The bounded O(σ_y/T) strain shift means finite-r values differ from the leading field by ~σ_y/σ_e (e.g. at r = 0.0613 mm, ~10%); the reported numbers are all leading (r → 0) coefficients as defined in the task.
- Screening endpoints would shift slightly if the checks were applied to the full finite-r field instead of the leading field (as stated, they are defined on the leading field).

Code: `solve.py` (calibration, Q2–Q4, full-law limit checks), `verify3.py` (displacement reconstruction, CTOD, strain checks, J contours).
