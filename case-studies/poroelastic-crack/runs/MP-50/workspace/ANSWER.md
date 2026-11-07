# A suddenly pressurised crack in a fluid-saturated porous solid

Reproduce with: **`python3 solve.py`** (prints the numbers written to `results.json`).

All results are for **(c)** ν = 0.2, ν_u = 0.4, α = 0.89, and are given *per unit
loading magnitude*, with **K in units of √(πL)** and **V in units of L²/G**
(the "G" is the shear modulus, so `V1_0 = 1.8849` means V = 1.8849 L²/G).

## Answers

| quantity | value | meaning |
|---|---|---|
| `K1_0`   | **1.0** | loading 1, t = 0⁺ |
| `K1_inf` | **1.0** | loading 1, t → ∞ |
| `V1_0`   | **1.8849555922** | loading 1, t = 0⁺ (undrained) |
| `V1_inf` | **2.5132741229** | loading 1, t → ∞ (drained) |
| `V1_short_coeff` | **1.8054066674** | A in V ≈ V1_0 + A·√(ct/L²) |
| `p_faces_0` | **0.7490636704** | pore pressure next to the faces, t = 0⁺ |
| `K2_0`   | **0.0** | loading 2, t = 0⁺ |
| `K2_inf` | **0.0** | loading 2, t → ∞ |
| `V2_inf` | **0.0** | loading 2, t → ∞ |

## Model and method

Infinite plane-strain Biot poroelastic body, slit crack −L ≤ x ≤ L on y = 0.
The two loadings are exactly the decomposition of the physical *suddenly
pressurised, permeable crack* (fluid pressure H(t) inside the crack):

* Loading 1 (σ_n = −1, p = 0) + Loading 2 (σ_n = 0, p = 1)  =  (σ_n = −1, p = 1),

and σ_n = σ_yy = −1 on a face is exactly the traction a fluid at pressure 1
exerts (traction on the solid face = −p_f·n, giving σ_n = −p_f). So **a negative
σ_n opens the crack**: the opening traction is τ = −σ_n = 1.

Three ingredients:

1. **Stress-intensity factor from the traction only** (Michel's theorem /
   Bueckner weight function). A cracked infinite body is *simply connected*;
   the divergence-free total-stress field with prescribed face tractions and
   zero field at infinity is *unique and independent of the elastic constants*,
   so the singular field is set by the crack-face tractions alone:
   K = τ√(πL) = −σ_n√(πL). This holds at t = 0⁺ (undrained) and t → ∞ (drained)
   alike, because both are stress fields with the *same* tractions; only the
   displacements (compliance) differ.
2. **Undrained / drained limits.** At t = 0⁺ the fluid has not moved, so the
   instantaneous response is elastic with **ν_u**; at t → ∞ the pore pressure
   has relaxed (loading 1: p ≡ 0), so the response is elastic with **ν**. The
   crack volume of a uniformly loaded elastic crack is V = π(1−ν_eff)τL²/G.
3. **Reciprocity (Betti) for the pore-pressure "eigenstrain"** gives the exact
   crack volume of the poroelastic solution and the short-time √t law.

## Derivation in brief

**Undrained pore pressure.** Tracing eq. (1),
2G·e_kk = (1−2ν)(σ_kk + 2αp). Setting ζ = 0 in eq. (2) (undrained) gives
σ_kk + α(1−2ν)p/(ν_u−ν) = 0, hence

    p_und = −(ν_u−ν)·σ_kk / [α(1−2ν)] .                              (U)

**Crack-face stress (loading 1).** For a uniformly loaded crack the stress on
the faces is hydrostatic: σ_xx = σ_yy = σ_n = −1 (standard pressurised-crack
solution: σ_xx+σ_yy = 4Re Φ, Φ(z) = (τ/2)(z/√(z²−L²) − 1), which on |x|<L gives
Re Φ = −τ/2, i.e. σ_xx = σ_yy = −τ). Hence **σ_kk = −2** on the faces. From (U):

    p_faces_0 = −(ν_u−ν)(−2)/[α(1−2ν)] = 2(ν_u−ν)/[α(1−2ν)] = 0.74906 .   (a)

**K.** Michel/weight function → K = −σ_n√(πL). Loading 1: σ_n = −1 ⇒
K1 = √(πL) at every time ⇒ **K1_0 = K1_inf = 1**. Loading 2: σ_n = 0 ⇒ **K2 = 0**
at every time. Consistency check: K1 + K2 = √(πL) is exactly the classical
hydraulic-fracturing result K = net-pressure·√(πL) for the physical crack.

**V for loading 1 (t = 0⁺ and t → ∞).** V = π(1−ν_eff)τL²/G with τ = 1:

    V1_0   = π(1−ν_u) = π·0.6 = 1.8849555922   (undrained, t = 0⁺),
    V1_inf = π(1−ν)   = π·0.8 = 2.5132741229   (drained, t → ∞).

**Short-time growth of V.** Write the eigenstrain ε*_ij = α(1−2ν)p/(2G)·δ_ij.
Betti's theorem for eigenstrain between the actual state A and the *drained
elastic reference crack* B (same face traction, p ≡ 0) gives, because the face
tractions coincide and the crack-face integrals reduce to the volumes,

    V^A − V^B = ∫ σ^B_kk(x,y) · α(1−2ν)p^A(x,y)/(2G)  dA ,          (R)

with σ^B_kk = 2(Re[z/√(z²−L²)] − 1) the trace of the unit-traction crack stress
(σ^B_kk = −2 on the faces). At small t the pore-pressure change from the
undrained initial value is confined to a boundary layer of thickness √(ct) at
the faces, following the 1-D diffusion profile
p^A − p_und = −p_und·erfc(y/2√(ct)) (both sides; p_und is *constant along each
face*, equal to (a)). Using σ^B_kk ≈ −2 in the thin layer,

    V(t) − V1_0 = 2·∫∫ (−2)·α(1−2ν)/(2G)·(−p_und) erfc(y/2√(ct)) dy dx
                = [α(1−2ν)p_und/G]·(2L)·2·(2√(ct)/√π)
                = 16(ν_u−ν)·L·√(ct)/(G√π)      (since α(1−2ν)p_und = 2(ν_u−ν)).

Hence, in units L²/G,

    V ≈ V1_0 + A√(ct/L²),   **A = 16(ν_u−ν)/√π = 16·0.2/√π = 1.8054066674** .  (b)

**Loading 2.** At t = 0⁺ the pore pressure is 1 only on the faces and 0 in the
interior (a measure-zero layer), so the solid is not yet loaded: K2_0 = 0,
V2_0 = 0. As t → ∞ the steady pore pressure is the harmonic field p = 1 on the
crack, p → 0 at infinity. A *harmonic* eigenstrain is compatible
(∇²p = 0 ⇒ no incompatibility), so it produces **no stress** in the
simply-connected body: the traction-free faces stay traction-free with zero
stress, hence **K2_inf = 0 and V2_inf = 0**. (Loading 2 is transiently non-zero
in between, but its two required limits are zero.)

## Checks made

**(1) Reciprocity coefficient.** Identity (R) must reproduce V1_0, which needs
∫(σ^B_kk)² dA = 2πL². Computed in elliptic coordinates
(z = L cosh(μ+iθ), dA = L²(sinh²μ+sin²θ)dμdθ, σ^B_kk = 2(coshμ sinhμ/D − 1),
D = sinh²μcos²θ + cosh²μsin²θ), Gauss–Legendre:

    [check] Int (sigma_kk)^2 dA = 6.283150   (target 2*pi = 6.283185)  ratio=0.999994

(refinement: 6.282867 → 6.283111 → 6.283150; larger `mumax` overflowed the
floats and was discarded). This moving number would be wrong for a wrong
coefficient in (R), a wrong branch of √(z²−L²), or a wrong Jacobian; the ratio
is 1 to 6 digits.

**(2) Short-time coefficient A.** Direct evaluation of (R) with the 1-D erfc
profile and the *exact* σ^B_kk(x,y) (not only its face value). As printed by
`solve.py`:

    [check] ct=1e-04  A_num=1.803599   A_formula=1.805407
    [check] ct=1e-05  A_num=1.803601   A_formula=1.805407
    [check] ct=1e-06  A_num=1.803602   A_formula=1.805407

A_num is independent of ct (confirming the √t law) and matches 16(ν_u−ν)/√π to
0.1%. On a refined grid the same integral reaches 1.805231 (0.01% agreement) and

    y=1.0e-06  int sigma_kk dx = -3.999600   (= -4: confirms sigma_kk = -2 on a face)
    y=1.0e-04  int sigma_kk dx = -3.999600
    y=1.0e-02  int sigma_kk dx = -3.999600

The residual is grid/tip truncation. This check would fail if the
boundary-layer profile, the factor 2 for the two faces, or the 1/(2G) in ε*
were wrong.

**(3) Physical consistency.** Loading1+Loading2 equals the permeable suddenly
pressurised crack, whose SIF must be the net traction times √(πL); here
K1+K2 = 1·√(πL) at both limits, matching the classical hydraulic-fracturing
result. A non-zero K2 would break superposition.

## Questions / assumptions (not fully pinned down)

* **Sign of "K".** I take σ_n = σ_yy = −1 (tension-positive σ_n) as the traction
  a fluid at pressure +1 applies to the face, so loading 1 *opens* the crack and
  K1 = +√(πL) (`K1_0 = K1_inf = +1`). Under the opposite convention the sign of
  K1 flips; magnitudes are unchanged. I report the opening sign.
* **Which "K".** I report the *total-stress* SIF (fixed by the traction boundary
  condition). An effective-stress SIF would include the bounded pore-pressure
  term; the weight-function result and the hydraulic-fracturing interpretation
  both use the total-stress SIF.
* **Face permeability.** Because *both* σ_n and p are prescribed on the faces for
  t > 0, the faces are permeable (drained in loading 1). This is what makes
  p_faces_0 a well-defined undrained-layer value and V2(t→∞) vanish.
* **`V2_0` and the loading-2 short-time law** were not requested; for
  completeness V2(0⁺) = 0 and V2 is transiently non-zero (returning to 0).
* **Diffusivity.** The short-time answer is expressed through the given
  c = 2κG(1−ν)(ν_u−ν)/[α²(1−2ν)²(1−ν_u)]; only the coefficient A of √(ct/L²)
  is reported, so κ enters only through c.