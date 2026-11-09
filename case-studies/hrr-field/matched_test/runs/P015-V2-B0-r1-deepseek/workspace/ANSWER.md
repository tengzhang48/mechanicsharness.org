# Answer: inferring T from C_zz, J

## Results (mm, N, MPa)

| Quantity | Value |
|---|---|
| **T** | **4816.47 MPa** |
| C_rr(60°) | **+3.9420×10⁻³ mm^1/2** |
| C_θθ(60°) | **+6.3318×10⁻⁴ mm^1/2** |
| C_γrθ(60°) = C_{2ε_rθ} | **+5.7311×10⁻³ mm^1/2** |
| C_y(90°) | **+9.1390×10⁻³ mm^1/2** |
| radius interval on θ=0 | **[0.01000, 0.061325] mm** (non-empty) |

Auxiliary quantities: K (tip amplitude) = 62.0740 MPa·mm^1/2; E_eff = T = 4816.47 MPa;
ν_eff = 0.486239; G_eff = 1620.356 MPa; κ_eff = 1.691358.

## Model chosen and why

**Field.** Mode I, plane stress (σ_zz = σ_xz = σ_yz = 0), traction-free faces, stationary straight
crack (β = π). As r → 0 the effective stress σ_e ∝ r^(−1/2) → ∞, so ⟨1 − σ_y/σ_e⟩ → 1 and the
given nonlinear total-strain law reduces to a **linear, isotropic** law. The leading stress field
is therefore the universal Williams mode-I field (stress exponent p = 1/2 fixed by J; its angular
functions are independent of the elastic moduli):

    σ_ij(r,θ) = (K/√(2πr)) f_ij(θ),
    f_rr = (5/4)cos(θ/2) − (1/4)cos(3θ/2),
    f_θθ = (3/4)cos(θ/2) + (1/4)cos(3θ/2),
    f_rθ = (1/4)sin(θ/2) + (1/4)sin(3θ/2).

**Constitutive reduction.** With the bracket → 1, the deviatoric compliance becomes
1/(2G_eff) = (1+ν)/E + (3/2)(1/T − 1/E), while the volumetric compliance (1−2ν)/E is unchanged.
This is an isotropic linear solid with the *same* bulk modulus and a modified shear modulus.
Its effective Young's modulus satisfies 1/E_eff = 1/(3G_eff) + 1/(9K), which here simplifies
**exactly to E_eff = T** (the 1/E terms cancel):

    1/E_eff = (2/3)[(1+ν)/E + 1.5(1/T−1/E)] + (1/3)(1−2ν)/E = 1/T.

**Treatment of ε_zz.** Plane stress leaves ε_zz free. In the reduced law (σ_zz = 0):

    ε_zz = −(ν/E)σ_kk + (3/2)(1/T−1/E)(−σ_kk/3) = −A_tr σ_kk,
    A_tr = ν/E + (1/2)(1/T−1/E) = (1/2)(1/T − 0.4/E).

On θ = 0, σ_kk = σ_xx+σ_yy = 2K/√(2πr), hence

    C_zz = lim √r ε_zz = −2A_tr K/√(2π) = −0.005000 mm^1/2.        (1)

**J.** J is path-independent for this (hyperelastic/deformation-theory) material; on a contour
shrinking to the tip the field is the above linear one, for which (verified numerically for
arbitrary isotropic moduli)

    J = K²/E_eff = K²/T = 0.800 N/mm.                            (2)

**Solving for T.** Eliminating K between (1) and (2):

    A_tr² T = π C_zz² / (2J) = 4.90874×10⁻⁵.

With A_tr = (1/2)(1/T − 0.4/E) this is a quadratic:
(1/T − 0.4/E)² T = 2π C_zz²/J. Its root with 0 < T < E is **T = 4816.47 MPa**, then
K = −C_zz√(2π)/(2A_tr) = 62.0740 MPa·mm^1/2, and J = K²/T = 0.8000 (consistent).

## Part 2 — strains at θ = 60°

Use the effective isotropic compliance (E_eff, ν_eff) with the angular stresses above. At 60°:
f_rr = 1.082532, f_θθ = 0.649519, f_rθ = 0.375000, and
C_rr = (K/√(2π))(f_rr − ν_eff f_θθ)/E_eff etc.:

    C_rr      = +3.9420×10⁻³
    C_θθ      = +6.3318×10⁻⁴
    C_γrθ     = (K/√(2π)) · 2(1+ν_eff) f_rθ / E_eff = +5.7311×10⁻³   (mm^1/2)

(Independent cross-check via the original three-term law with the same f's gives identical values.)

## Part 3 — C_y(90°)

The mode-I plane-stress displacement field is u_y = (K/2G_eff)√(r/2π) sin(θ/2)(κ_eff − cosθ)
with κ_eff = (3−ν_eff)/(1+ν_eff) = 1.691358. At θ = 90°:

    C_y(90°) = lim u_y/√r = K κ_eff / (4 G_eff √π) = +9.1390×10⁻³ mm^1/2.

This was reproduced to 12 digits by integrating the strains (u_r = 2C_rr r^1/2, then
U_θ′ = C_θθ_ang − U_r), a check independent of the measured C_zz.

## Part 4 — radius screen on θ = 0

**σ_e ≥ 10σ_y = 100 MPa.** On θ=0, σ_xx = σ_yy = K/√(2πr), σ_zz = 0, and σ_e = K/√(2πr).
Thus r ≤ K²/(2π·100²) = **0.0613252 mm**.

**Largest principal strain magnitude ≤ 0.0500.** At θ=0 the principal strains are the double
in-plane strain and ε_zz. With f_rr = f_θθ = 1, the in-plane value is
(K/√(2π))(1−ν_eff)/E_eff · r^(−1/2) = 0.0026408/√r, and from (1)
ε_zz = C_zz/√r = −0.005/√r. So |ε_zz| ≤ 0.05 for r ≥ 0.005/0.05 = **0.010000 mm** (the binding
strain; the in-plane one is 0.0264 there). Hence r ≥ 0.01000 mm.

**r ≤ 0.200 mm** is not binding (0.0613 < 0.200).

Combining: the interval is **[0.01000, 0.061325] mm**, non-empty (length 0.051325 mm).

## What the screens establish; crack initiation

The three screens bound the annulus where the leading field is internally consistent: σ_e ≥ 10σ_y
keeps the dropped ⟨1−σ_y/σ_e⟩ factor within 10% (so the linear reduction is accurate), the strain
cap keeps infinitesimal kinematics valid, and r ≤ 0.200 mm keeps outer-body corrections negligible.
They are defined consistency checks on the asymptotic field, **not** error bars for the full
finite-r solution.

The supplied data do **not** determine crack initiation. J = 0.800 N/mm and the fields give the
tip loading and K, but a propagation/initiation threshold (a toughness J_c, K_IC or notch criterion)
is not provided; J alone is a loading parameter, not a failure criterion. The σ_y and T values set
the constitutive response only.

## Verification (beyond substitution into C_zz)

1. **J from the field.** Numerically evaluating J = ∫_Γ(W N_x − t_i u_i,x) ds on circles of radius
   0.001–0.05 mm with E_eff, ν_eff, K (Cartesian components, incl. polar-basis derivatives) gives
   J = 0.8000 for all radii — confirms J = K²/E_eff and the plane-stress convention.
2. **Independent displacement integration** reproduced C_y(90°) to 12 significant digits.
3. **Provided checks.** `knowledge/checks/crack_tip_checks.py` run on the Williams field with
   ν_eff: traction-free faces, mode-I symmetry, angular and Cartesian equilibrium, Saint-Venant
   compatibility, finite energy, p = 1/2, crack-opening and strain→displacement compatibility all
   PASS.
4. **Algebraic identity** E_eff = T and the closed-form quadratic root were cross-checked against
   the numerical solve.

## Missing or uncertain

- Path independence of J (equivalence of the externally supplied J and the near-tip J) is assumed;
  it holds for a hyperelastic/deformation-theory material, which the given law is taken to be.
- The reduction drops O(σ_y/σ_e) corrections; this is asymptotically exact as r → 0 and the screen
  σ_e ≥ 10σ_y quantifies where it is good to ~10%.
- C_zz is supplied to 5 significant figures, which dominates the numerical uncertainty (~0.1% in T).
- Plane-stress idealization is taken from the stated 2-D constraint (σ_zz = 0, ε_zz free).
