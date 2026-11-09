# Answer

## 1. Results

| # | Quantity | Value |
|---|---|---|
| 1 | **T** | **4816.47 MPa** (solve: T = 4816.4722; ν* = 0.4862387, K = 62.0740 MPa·mm^1/2) |
| 2 | θ = 60°: C_rr | **0.0039420 mm^1/2** |
| 2 | θ = 60°: C_θθ | **0.00063318 mm^1/2** |
| 2 | θ = 60°: C_γrθ (= lim √r·2ε_rθ) | **0.0057311 mm^1/2** |
| 3 | C_y(90°) = lim u_y(r,90°)/√r | **0.0091390 mm^1/2** |
| 4 | screening interval on θ = 0 | **r ∈ [0.01000, 0.06133] mm** (non-empty) |

Exact closed forms used (so the numbers are reproducible): with ν* = 1/2 + (ν − 1/2)T/E,
K = √(JT), and f = K/(T√(2π)):

- C_zz(0) = −2ν*√(J/T)/√(2π) = −0.005  →  T is the unique root in (0, E) (LHS strictly
  decreasing; equivalently the quadratic (1/2 − aT)² = c²T, a = (1/2 − ν)/E,
  c = 0.005√(2π)/(2√J)).
- At 60° (where σ_rr, σ_θθ, σ_rθ have angular values 5√3/8, 3√3/8, 3/8):
  C_rr = (√3/8)(5 − 3ν*)f, C_θθ = (√3/8)(3 − 5ν*)f, C_γrθ = (3/4)(1 + ν*)f.
- C_y(90°) = (3 − ν*)K/(2√π T) = 2·C_rr(90°) (two independent routes agree to 10 digits).
- Screens on θ = 0 (leading field): σ_e = K/√(2πr) ≥ 10σ_y gives
  r ≤ JT/(2π(10σ_y)²) = 0.061325 mm; max |principal strain| = |ε_zz| = 2ν*σ/T = 0.005/√r ≤ 0.05
  gives r ≥ (0.005/0.05)² = 0.010000 mm; cap r ≤ 0.200 mm is inactive. Endpoints: at
  r = 0.01, σ_e = 247.6 MPa; at r = 0.06133, |ε_zz| = 0.0202.

## 2. Model chosen and why

Loading is monotonic and proportional at the material points, so total (deformation) theory is
the correct constitutive frame and a strain-energy/J exists (knowledge TOPIC §2, §6 item 6).
As r → 0 the stress is singular, σ_e → ∞, so the Macaulay factor ⟨1 − σ_y/σ_e⟩ → 1 and the law
reduces to a **linear isotropic compliance**:

ε_ij ≈ [(1+ν)/E + (3/2)(1/T − 1/E)]σ_ij + [−ν/E − (1/2)(1/T − 1/E)]σ_kk δ_ij,

which is Hooke's law with effective constants **E\* = T** and **ν\* = 1/2 + (ν − 1/2)T/E**
(check: uniaxially ε = σ/T + σ_y(1/E − 1/T), a bilinear curve with tangent modulus T; the
plastic term is incompressible, hence ν* → 1/2 as T → 0). The deviations from this linear law
are O(σ_y/σ_e) = O(√r), i.e. bounded strains that do not affect any √r-coefficient.

The leading field is therefore the mode I Williams field for a crack (β = π): exponent
p = 1/2 (fixed by finite, nonzero J — the linear case of the HRR argument; energy admissibility
p < 1 satisfied), with angular functions independent of ν*, under plane stress (σ_zz = σ_xz =
σ_yz = 0). This is exactly Hutchinson's piecewise-linear result (entry, eq. 42): amplitude
factor √(E_tan/E), i.e. here K_eff = σ∞√(πT) from J = πσ∞² — equivalently **J = K²/E\* = K²/T**
(plane stress), so **K = √(JT)**.

**Treatment of ε_zz**: the constraint is plane stress, so σ_zz = 0 and ε_zz is the free
kinematic response: ε_zz = −ν*σ_kk/E* = −ν*(σ_xx + σ_yy)/T to leading order. On θ = 0,
σ_kk = 2K/√(2πr), giving the calibration equation above. ε_zz carries the through-thickness
Poisson response of the combined elastic + (incompressible) plastic compliances.

## 3. Assumptions changed

None of the stated assumptions was changed. Two judgement calls, both stated:
(a) the *leading* field is treated as exactly the effective-linear field — this is exact for
the r → 0 coefficients by definition; applying it unmodified out to r = 0.200 mm relies on the
problem's statement that matching corrections are negligible there (it overrides Hutchinson's
own reservation that a single-term bilinear field is "suspicious"); (b) deformation theory is
used as given (proportional loading is assumed), not flow theory.

## 4. Verification (solve.py, all run)

1. **Independent J check (does not use C_zz):** reconstructed the full leading field in real
   units (K-field stresses, strains via Hooke(T, ν*), displacements via the plane-stress mode I
   formulas) and evaluated J = ∫(W N_x − t_i u_{i,x})ds on circular arcs with u_{i,x} by finite
   differences in Cartesian components: **J = 0.800000 N/mm at R = 0.02, 0.05, 0.1, 0.2 mm**
   — reproduces the supplied J from the amplitude K = √(JT), validating the angular functions,
   the displacement field and the J–K relation together.
2. **Full-law convergence:** √r·ε_zz(r, 0) from the complete nonlinear law approaches −0.005000
   monotonically (−0.004727 at r = 0.02; −0.004939 at r = 0.001) with error = 0.00193√r and
   extrapolated r → 0 intercept −1.8e−19; likewise the 60° coefficients converge to the closed
   forms (e.g. C_rr: 0.003913 at r = 5e−4 vs 0.003942). Confirms the constitutive reduction
   (Macaulay factor → 1) and that the reported numbers are the r → 0 limits.
3. **Knowledge checks** (`knowledge/checks/crack_tip_checks.py`, law = "linear", ν = ν*,
   plane stress, p = 0.5, mode I): all pass — faces traction-free, mode I symmetry, angular
   and Cartesian equilibrium, Saint-Venant compatibility, finite energy, crack exponent 1/2,
   faces open, displacement compatibility, and the field's own J (6.28319 = 2π in the E\* = 1
   normalization, matching the input K-normalization to 7e−14).
4. C_y(90°) obtained by two independent routes (displacement formula vs u_y = u_r = 2C_rr√r
   integration along the ray): identical to 10 digits. T obtained by two independent solves
   (Brentq and the closed-form quadratic): identical.

## 5. What the screens establish; what is missing or uncertain

The screens delimit where the **leading field is a self-consistent description**: σ_e ≥ 10σ_y
puts the point deep in the plastic range where dropping σ_y/σ_e in the reduction costs ≲ 10%
(relative error of the Macaulay factor) and deformation ≈ flow theory under proportionality;
max |principal strain| ≤ 0.05 keeps infinitesimal kinematics defensible (~5% strain); r ≤ 0.2 mm
stays inside the annulus where matching to the surrounding body is given to be negligible. So
[0.01, 0.0613] mm is the computed window in which this asymptotic field can be quoted with its
own stated accuracy. They are defined checks on the leading field, **not** error estimates of
the full finite-r solution, and they say nothing about failure.

**Crack initiation is not determined.** The supplied data (constitutive constants, J, and the
field) are field/amplitude information; initiation requires a fracture criterion (critical
K, J, CTOD, or a damage model), which is not part of the problem. Per the knowledge topic
(§6 item 9), J is a field parameter, and energy release becomes a fracture parameter only
through an additional toughness assumption.

Uncertainties: (i) the leading field breaks down near the crack faces, where σ_e < σ_y
(plasticity switches off; elastic sector/transition at angles where cos(θ/2) ~ σ_y√(2πr)/K);
(ii) actual finite-r values differ from the leading coefficients by O(σ_y/σ_e) ≈ O(√r) —
about 2% at r = 0.01 mm and 18% at r = 0.2 mm for θ = 0 — acceptable here only because the
problem guarantees negligible corrections for r ≤ 0.2 mm; (iii) real 3D fronts, thickness, and
possible face contact are outside the 2D idealization (§6 items 5, 8).
