# ANSWER — Infer T from a transverse strain field

## 1. Requested numbers

| Quantity | Value |
|---|---|
| **T** | **4816.47 MPa** (≈ 4.816 GPa; 1% band 4768–4865) |
| C_rr(60°) = lim √r ε_rr | **3.9420×10⁻³ mm^1/2** |
| C_θθ(60°) = lim √r ε_θθ | **6.3318×10⁻⁴ mm^1/2** |
| C_γrθ(60°) = lim √r (2ε_rθ) | **5.7311×10⁻³ mm^1/2** |
| C_y(90°) = lim u_y/√r | **9.1390×10⁻³ mm^1/2** |
| Screening interval on θ = 0 | **r ∈ [0.0100, 0.0613] mm** (non-empty; endpoints 0.010000 and 0.061325 mm) |

Auxiliary: ν_eq = 0.486239, K_I = √(J·T) = 62.074 MPa·mm^1/2, K_I/√(2π) = 24.7639 MPa.

## 2. Model chosen and why

**Total-strain (deformation) theory, leading separable term.** Loading is quasistatic, monotonic and
proportional at the material points, so deformation theory is exact for this material and W = ∫σ:dε is a
potential; the crack is stationary with traction-free faces parallel to x, the membrane homogeneous, so
**J is path-independent**.

**Constitutive reduction (including ε_zz).** With σ ~ r^(−p) and p = 1/2 (below), σ_e = √(3/2 s:s) ~ r^(−1/2) → ∞
as r → 0, so ⟨1 − σ_y/σ_e⟩ → 1 and the term −(σ_y/σ_e)s_ij is O(r⁰) in strain — subleading. The leading-order
law is then exactly linear and isotropic:

  ε_ij = (1+ν)/E σ_ij − ν/E σ_kk δ_ij + (3/2)(1/T − 1/E) s_ij  ⟺  Hooke with bulk modulus K = E/(3(1−2ν))
  (unchanged) and deviatoric coefficient A = 3/(2T) + (ν − ½)/E, i.e. an **equivalent isotropic material**

  **E_eq = T (exact identity 1/E_eq = 1/(9K) + (2/3)A = 1/T),  ν_eq = ½ − (1−2ν)T/(2E) = 0.486239.**

ε_zz is unconstrained (σ_zz = 0, plane stress) and follows from this law: ε_zz = −(ν_eq/T)(σ_xx+σ_yy); it does
not feed back on the 2D field. The measured C_zz is precisely this follower strain, which is what makes T observable.

**Field.** Finite energy + finite non-zero J for a crack ⇒ p = 1/2 (Williams; HRR n→1 limit). The leading stress
field is the standard mode-I Williams field σ_ij = (K_I/√(2πr))·f_ij(θ) (angular functions independent of the
elastic constants); strains follow from the equivalent linear law. Symmetric tension ⇒ pure mode I (K_II = 0).
Amplitude from J evaluated on a vanishing circle, where the leading field dominates and subleading terms drop out:
**J = K_I²/E_eq = K_I²/T** (plane stress), consistent with Hutchinson 1968 eq. 42 for the piecewise-linear
material (tip-amplitude factor √(E_tan/E) relative to the remote elastic K).

**Determining T.** On θ = 0, σ_xx = σ_yy = K_I/√(2πr), so
C_zz = −2ν_eq K_I/(√(2π) T) = −√(2/π)·√J·ν_eq/√T with ν_eq = ½ − (1−2ν)T/(2E).
Solving |C_zz| = 0.00500 (brentq, `solve.py`) gives **T = 4816.47 MPa**.

Then, with fac = K_I/(√(2π)T): C_rr = (f_rr − ν_eq f_θθ)fac, C_θθ = (f_θθ − ν_eq f_rr)fac,
C_γrθ = 2(1+ν_eq) f_rθ fac, where at θ = 60° the polar angular functions are f_rr = 5√3/8 = 1.082532,
f_θθ = 3√3/8 = 0.649519, f_rθ = 3/8 = 0.375 (literature values, reproduced exactly).
Displacements (no rigid translation; u ~ √r → 0 at the tip): u = (K_I(1+ν_eq)/T)√(r/2π)·{c(θ/2)[κ−1+2sin²(θ/2)],
s(θ/2)[κ+1−2cos²(θ/2)]}, κ = (3−ν_eq)/(1+ν_eq) = 1.691358; C_y(90°) = K_I(1+ν_eq)κ/(2√π T) = 9.1390×10⁻³.

## 3. Assumptions changed

None of the problem's assumptions were changed. Interpretation choices I made: (a) pure mode I from the stated
symmetric tension; (b) J and C_zz treated as exact inputs, and the "negligible matching corrections for
r ≤ 0.200 mm" taken as given (so only the leading term is needed there); (c) displacements taken as the standard
K-field particular solution (no rigid translation); (d) the formal r → 0 coefficient used as defined, even though
the screens below restrict where the field is physically usable.

## 4. Verification (beyond substituting into C_zz)

1. **Independent check suite** (`knowledge/checks/crack_tip_checks.py`, run on my angular functions with
   law = "linear", plane = "stress", ν = ν_eq): traction-free faces 3×10⁻²⁴; mode-I symmetry 10⁻¹⁶; equilibrium
   (angular ODEs) 10⁻¹⁰ and (Cartesian FD) 10⁻⁸; Saint-Venant compatibility 4×10⁻⁷; exponent p = 1/2 exactly; faces
   open; displacement compatibility 3×10⁻¹²; and the field's own J (computed by the checker's Cartesian machinery)
   = 3853.18 = J·T = K_I² to 10⁻¹³ — an independent confirmation of J = K_I²/T with my normalization.
2. **J with the true nonlinear W** (deformation-theory potential including the σ_y cutoff), integrated on circles:
   0.800000 N/mm at r = 10⁻⁵–10⁻³ mm; 0.7999 at 0.0613; 0.7995 at 0.2 — i.e. subleading corrections to the
   leading-field J are < 0.06 % throughout the screening annulus.
3. **C_y three independent ways**: integrating the polar strains (residual 8×10⁻¹³), the closed form, and the line
   integral u_y(0,y) = 2C_εyy(90°)√y along x = 0 — all 9.139019×10⁻³.
4. Closed-form displacements reproduce the law strains by finite differences to ~10⁻¹²; full nonlinear law gives
   √r·ε_zz(r,0) → −0.00500 (−0.00481 at r = 0.01, the subleading bracket term decaying like √r).

## 5. Screening interval (θ = 0), and what the screens establish

- σ_e = K_I/√(2πr) ≥ 10σ_y = 100 MPa ⇒ r ≤ K_I²/(2π(10σ_y)²) = 0.061325 mm.
- Principal strains on θ = 0 are (ε, ε, ε_zz) with ε = 2.6415×10⁻³/√r and ε_zz = −5.000×10⁻³/√r (largest);
  ≤ 0.0500 ⇒ √r ≥ 0.005/0.05 ⇒ r ≥ 0.010000 mm.
- r ≤ 0.200 mm (given matching) — not binding.

**Interval [0.0100, 0.0613] mm.** The σ_e screen keeps the material deeply plastic (⟨1−σ_y/σ_e⟩ ≥ 0.90), which is
the self-consistency condition for the linear-hardening leading field; the strain screen keeps infinitesimal
kinematics tolerable (≤ 5%); r ≤ 0.2 is where matching corrections are stated negligible. Together they bracket
the annulus where the leading field is a usable approximation of the full solution — defined checks, not error
estimates of the finite-r solution.

**Do the data determine crack initiation?** No. J = 0.800 N/mm fixes the field amplitude only. Initiation needs a
fracture criterion — a critical K/J (toughness), a critical stress/strain, or a damage model — none is supplied
(the material model has no failure surface). Also absent: thickness (plane-stress 2D needs the field zone ≫
thickness), large-strain/blunting physics near the tip (process-zone scale ~J/σ_y ≈ 0.08 mm, comparable to the
annulus), and any statement about load history beyond the current J.

## 6. Missing / uncertain

- Only the leading term is computed; higher-order terms (elastic-plastic transition field, T-stress-like terms)
  are unquantified outside the annulus but stated negligible for r ≤ 0.2 mm.
- If the real constraint were not plane stress, or loading not exactly proportional, the reduction J = K_I²/T and
  ν_eq would change (plane strain would give J = (1−ν_eq²)K_I²/T and a different C_zz equation).
- The 2% tolerances are easily met by the formulation itself; the dominant uncertainty is model validity, not
  arithmetic. Files: `solve.py` (full solution + checks), `run_output.txt` (log), `PROCESS.md`.
