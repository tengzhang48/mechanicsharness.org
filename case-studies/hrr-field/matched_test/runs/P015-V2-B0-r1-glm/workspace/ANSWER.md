# ANSWER — Infer a material stiffness from a transverse strain field

Units: mm, N, MPa. Script: `solve.py` (run `python3 solve.py`); log: `PROCESS.md`.

## 1. Requested numbers

| Quantity | Value |
|---|---|
| **T** | **4816.47 MPa** (≈ 4.816 GPa) |
| C_rr(60°) = lim √r ε_rr | **+0.003942 mm^1/2** |
| C_θθ(60°) = lim √r ε_θθ | **+0.000633 mm^1/2** |
| C_γrθ(60°) = lim √r (2ε_rθ) | **+0.005731 mm^1/2** |
| C_y(90°) = lim u_y(r,90°)/√r | **+0.009139 mm^1/2** |
| Screen interval on θ = 0 | **r ∈ [0.02336, 0.06133] mm** (non-empty) |

Derived quantities: effective leading-field moduli E* = T = 4816.47 MPa, ν* = 1/2 − (1/2 − ν)T/E = 0.48624; stress intensity K* = √(J T) = 62.074 MPa·mm^1/2.

## 2. Model chosen and why

**Constitutive reduction.** Near the tip σ_e = K*/√(2πr)·(shape) → ∞ as r → 0, so ⟨1 − σ_y/σ_e⟩ = 1 − σ_y/σ_e → 1 and the omitted part is bounded (O(σ_y), i.e. relative size σ_y/σ_e = O(√r)). The total-strain law reduces asymptotically to a **linear** isotropic law

ε_ij = a σ_ij + b σ_kk δ_ij,  a = (3/2)/T + (ν − 1/2)/E,  b = (1/2 − ν)/E − 1/(2T),

which is exactly Hooke's law with **E\* = T and ν\* = 1/2 − (1/2 − ν)T/E = 0.4862** (verified algebraically in the script to machine precision: a = (1+ν*)/T, −ν*/T = b). Physically: beyond yield the uniaxial tangent is T, and the plastic deviatoric term drives the response toward incompressibility (ν* → 1/2 as T/E → 0).

**Field.** With the asymptotically linear law the tip eigenproblem is Williams's; J is path-independent (a strain-energy density exists for this J2 total-strain — Hencky — material, and proportional loading is given, so deformation theory is exact here), the flanks are traction-free and parallel to x, so finite nonzero J fixes the exponent p = 1/2 (checked in §4). The loading is symmetric → pure mode I. Hence the leading field is the **plane-stress mode I K-field with moduli (E*, ν*) = (T, ν\*)**:

- σ_ij = (K*/√(2πr)) f_ij(θ), σ_zz = 0;
- ε_ij from the compliance (a, b) above — including **ε_zz = b(σ_xx+σ_yy)**, which is *free* under the σ_zz = 0 constraint and carries the measured C_zz (note b < 0 → ε_zz < 0 ahead of the tip);
- u_i = (1+ν*)K*/T · √(r/2π) × standard mode I angular functions with κ = (3−ν*)/(1+ν*) = 1.6914, no rigid translation.

**Amplitude from J.** On a vanishing circular arc the field is the asymptotic linear one, so J = K*²/E* with E* = T (plane stress; the (1−ν²) factor belongs to plane strain only): **K* = √(J T)**. This plus the measured C_zz(0) = 2b·2K*/√(2π) = −0.00500 determines T: the left side is strictly increasing in T on (0, E) (derivative ∝ 1/(2T) + (1/2−ν)/E > 0, from −0.357 at T→1 to −8.1e−4 at T→E), so the root T = 4816.47 MPa is unique. Then C_rr, C_θθ, C_γrθ at 60° are the polar components of the strain coefficients (rotation with e_r = (cosθ, sinθ), e_θ = (−sinθ, cosθ); C_γrθ is the *engineering* shear coefficient as requested), and C_y(90°) = (1+ν*)/T · K*/√(2π) · sin45°·[κ+1−2cos²45°] = (1+ν*)K*κ/(2T√π).

**Screening interval (θ = 0).** There σ_xx = σ_yy = p = K*/√(2πr), σ_xy = 0, so σ_e = p exactly (deviator (p/3, p/3, −2p/3)); principal strains are ε_xx = ε_yy = a·p = 3.0857e−4·p (largest in magnitude; |ε_zz| = 2.019e−4·p). Screens: σ_e ≥ 100 MPa ⇒ r ≤ K*²/(2π·100²) = 0.061325 mm; a·p ≤ 0.05 ⇒ r ≥ (a K*/(√(2π)·0.05))² = 0.023357 mm; r ≤ 0.200 mm (not active). Interval **[0.02336, 0.06133] mm**.

## 3. Assumption changes

None of the stated problem assumptions was changed. Two reductions were *added* and justified rather than assumed: (i) ⟨1−σ_y/σ_e⟩ → 1 as r → 0, making the law asymptotically linear (exact-law numerics below confirm convergence); (ii) J evaluated on the asymptotic field, valid because J is path-independent for this Hencky material and bounded corrections vanish on shrinking contours. I also assumed pure mode I (from the stated symmetric loading) and no mode II component.

## 4. Verification (beyond substituting into C_zz)

1. **Exact-law limit (not the asymptotic reduction):** evaluating the *full* nonlinear law at finite r on θ = 0 gives √r·ε_zz = −0.004807 (r = 10⁻²), −0.004939 (10⁻³), −0.004981 (10⁻⁴), −0.004998 (10⁻⁶), −0.005000 (10⁻⁸) → the measured −0.005000. Same check at θ = 60° reproduces C_rr, C_θθ, C_γrθ to 6 digits at r = 10⁻⁸.
2. **Provided checks** (`knowledge/checks/crack_tip_checks.py`, law="linear", nu=ν*, plane="stress", mode I, β=π, p=1/2, amplitude √J with E*=1): **all pass** — traction-free faces, mode I symmetry, angular equilibrium (8e−10), Cartesian equilibrium by finite differences (1e−8), Saint-Venant compatibility (4e−7), finite energy, exponent p = 1/2 fixed by J, faces open, displacement–strain shear relation (2e−12), and the field's own J = 0.800 at r = 0.5, 1, 2 (matches the supplied J to 1e−13). These checks also caught (and I fixed) a sign error in my first polar rotation of the stress shapes.
3. **Displacement closure:** differentiating the analytic u field numerically reproduces the compliance strains to 7.5e−13 at (r, θ) = (1, 60°); and C_y(90°) = 2·(ε_rr coefficient at 90°), since u_y(90°) = u_r(90°) = 2√r·Ê_rr — the two independent routes agree (0.009139). This cross-check caught a wrong prefactor K/2E vs (1+ν)K/E in my first u_y formula.

## 5. What the screens establish; crack initiation

The three screens define the annulus in which the leading field is a *self-consistent* description: σ_e ≥ 10σ_y keeps the material deeply plastic (the neglected factor σ_y/σ_e ≤ 0.1, so the asymptotic linear reduction is uniformly accurate to ≲10% and far from the elastic–plastic boundary); max principal strain ≤ 0.05 keeps infinitesimal kinematics tenable (5% is already marginal); r ≤ 0.200 mm is where matching to the surrounding body is stated negligible. They are defined checks on the leading field, not error estimates of the full finite-r solution; outside the annulus (r < 0.0234 mm: kinematics; r > 0.0613 mm: yield-screen/log correction ~10%) the single leading term alone should not be used.

**Crack initiation is not determined.** The constitutive law plus J fix the *field* (an amplitude parameter), not a failure condition: no toughness (J_IC), critical CTOD, or critical-strain/damage criterion is supplied, and damage is explicitly excluded. Predicting initiation requires an additional fracture criterion; the topic file (§6, item 9) makes the same point.

## 6. Missing / uncertain

- The supplied C_zz(0) and J are taken exact as stated; T inherits their precision (the 1% tolerance is generous — a 1% change in C_zz moves T by ≈1%).
- Uniqueness of T relies on monotonicity shown in §2; verified numerically over (0, E).
- J-dominance (one-parameter characterization) is assumed per the problem's statement that matching corrections are negligible for r ≤ 0.200 mm; the plastic-zone size implied by the leading field (σ_e = σ_y near r ≈ 6 mm) is far outside that window, so physically J-dominance of the outer problem must come from the external calculation, not from this field.
- Plane stress is a stated 2D constraint, not a 3D claim; a real membrane's tip field needs the plane-stress zone to exceed the thickness (topic §6, item 8).
