# Mode III crack-tip field in a power-law shear material — results

Units: mm, N, MPa. Material: γ = (γ₀/τ₀)(τ/τ₀)^(n−1) τ with τ₀ = 80 MPa, γ₀ = 0.002, n = 3.5
(pure power law, no elastic term). J = 0.120 N/mm. Solver: `mode3.py` (IVP shooting + quadrature).

## 1. The seven requested numbers

| Quantity | Value |
|---|---|
| Δw(r = 0.0200 mm) | **1.3181 × 10⁻³ mm** |
| Δw(r = 0.0800 mm) | **1.7937 × 10⁻³ mm** (ratio (4)^{2/9} = 1.36079 ✓) |
| τ_rz(r = 0.0200 mm, θ = 60°) | **+38.59 MPa** |
| τ_θz(r = 0.0200 mm, θ = 60°) | **+138.0 MPa** |
| r: \|γ\| = 0.0100 at θ = 0° | **0.04689 mm** |
| r: \|γ\| = 0.0100 at θ = 90° | **0.02507 mm** |
| r: \|γ\| = 0.0100 at θ = 180° (face limit) | **0.01340 mm** |

Signs: τ_θz(θ = 0) > 0 (τ_yz > 0 ahead of the tip) and w(r, +π) > 0, consistent with the upper
face moving in +z; Δw > 0. All intersection radii are ≪ 0.30 mm, so the leading field applies.

## 2. Model chosen and why

Governing model: the **separable anti-plane (mode III) asymptotic field** for a traction-free
crack (β = π) in the pure power-law shear law. With τ = r^(−p) A_T T(θ) = r^(−p) A_T (T_r, T_θ)
and w = r^s F(θ):

- equilibrium: (1 − p) T_r + T_θ′ = 0;
- compatibility of γ = ∇w with γ = (γ₀/τ₀ⁿ) |τ|^(n−1) τ: (m^(n−1) T_r)′ = s m^(n−1) T_θ,
  m = |T|, s = 1 − np (so γ ∝ r^(−np), w ∝ r^s);
- BCs from the stated symmetry and loading: T_r(0) = 0 (w odd in y ⇒ τ_xz = 0 on θ = 0) and
  T_θ(π) = 0 (traction-free faces); T_θ(−θ) = T_θ(θ), T_r odd.

Why this model: the field is the leading term of the intermediate asymptotics; the flanks are
parallel to x and traction-free, no body forces, a strain energy exists (W = n/(n+1) τ·γ for the
power law), so **J is path-independent and fixes the exponent**: a separable field gives
J ∝ r^(1−p(n+1)), and finite nonzero J requires **p = 1/(n+1) = 2/9** (TOPIC.md §3 argument,
applied here in anti-plane form). The nonlinear two-point problem was solved by shooting on p
with T_r(0) = 0, T_θ(0) = 1 (`solve_ivp`, rtol 1e-12): the shooting eigenvalue is
**p\* = 0.222222222222 = 1/(n+1) to 2.8 × 10⁻¹⁷** — the eigenvalue condition and the J condition
agree, which is a strong internal consistency check.

How loading fixes the magnitude: with the angular functions normalized (T_θ(0) = 1), the field's
own energy flux through a circular contour is J = (γ₀/τ₀ⁿ) A_T^(n+1) I_n with
I_n = 2∫₀^π [(n/(n+1)) m^(n+1) cosθ − m^(n−1) T_r² cosθ + m^(n−1) T_r T_θ sinθ] dθ = 2.01960.
Setting J = 0.120 N/mm gives the amplitude **A_T = 64.193 MPa·mm^(2/9)** (τ = A_T r^(−p) T,
w = r^s F, F = A_Tⁿ (γ₀/τ₀ⁿ) m^(n−1) T_r / s). Then Δw = 2 r^s F(π), stresses from
τ_rz = A_T r^(−p) T_r, τ_θz = A_T r^(−p) T_θ, and \|γ\| = r^(−np) (γ₀/τ₀ⁿ)(A_T m)^n.

## 3. Assumptions changed

None of the problem's physical assumptions were changed. One modelling note: the supplied
`knowledge/checks/crack_tip_checks.py` covers in-plane (mode I/II) fields only, so it does not
apply to anti-plane shear; I wrote equivalent mode III checks instead (below). Deformation theory
is legitimate here because, inside the dominance zone, monotonic J makes the stress at a fixed
material point proportional (TOPIC.md §6 item 6), as the problem's "quasistatic, monotonic, fixed
loading pattern" states.

## 4. Verification

1. **n = 1 regression against the exact mode III elastic field** (whole pipeline, same code):
   p\* = 1/2 (err 1e-16), T = (sin θ/2, cos θ/2) to 1e-15, I₁ = π to 4e-15, and Δw(r = 0.05)
   matches the closed form w = (K_III/μ)√(2r/π) sin(θ/2) with K_III = √(2μJ), μ = τ₀/γ₀, to
   1 × 10⁻¹⁵ relative. This validates the ODEs, the J integral and its sign conventions end to end.
2. **Boundary residual**: \|T_θ(π)\| = 1.3 × 10⁻¹⁶ (traction-free faces), and F odd in θ as required.
3. **Energy flux**: J recomputed from the assembled n = 3.5 field in **Cartesian components**
   (t_z from Cartesian stresses, w_,x by finite differences in x, W = n/(n+1) τ·γ) on circular
   contours R = 0.005, 0.02, 0.05, 0.15 mm: J = 0.1200000 N/mm on every one (target 0.120).
4. **Equilibrium** in Cartesian components by finite differences at sample points: relative
   residual ≤ 7 × 10⁻⁶ (finite-difference limited).
5. Eigenvalue from shooting equals the J-fixed value 2/9 to 3 × 10⁻¹⁷ (item 2 above), and
   Δw(0.08)/Δw(0.02) = 4^(2/9) = 1.36079 exactly, as the separable form requires.

## 5. Missing / uncertain

- Only the **leading term** is used; valid for r ≪ 0.30 mm (all requested radii are ≤ 0.08 mm;
  the problem states geometry beyond 0.30 mm is a negligible correction). Higher-order terms are
  not computed.
- Small strain fails as r → 0 (\|γ\| → ∞), so the numbers lose physical meaning at sufficiently
  small r; the requested points have γ = O(10⁻²), still small.
- The answer to the final question: **No** — these results (even with J fully fixed) characterize
  the crack-tip driving force only, and the load at which the slit starts to advance additionally
  requires a material failure criterion (e.g. a critical J, CTOD, or damage-based condition) that
  the problem does not supply.
