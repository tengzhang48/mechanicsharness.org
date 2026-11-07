# ANSWER — a suddenly pressurized crack in a fluid-saturated porous solid

**Reproduce:** `python3 solve.py` (prints the values below, rewrites `results.json`,
then runs independent numerical checks; the run used here is pasted in `run_output.txt`).

## Model and method, and why

The solid is a Biot poroelastic plane with a Griffith crack of half-length L. Two facts
make the problem solvable in closed form at the times requested:

1. **Instantaneous response = undrained elasticity.** A suddenly applied load leaves no
   time for fluid content to change (ζ = 0 at t = 0⁺), so the solid responds as an
   elastic solid with the *undrained* Poisson ratio ν_u and the same G. At t → ∞ the pore
   pressure has equilibrated; for loading 1, p → 0 everywhere (faces drained, no source),
   so the response is the *drained* elastic one (G, ν).

2. **Griffith crack formulas** (unit opening traction σ₀, operative Poisson ratio ν′):
   - K_I = σ₀√(πL) — set by the traction boundary condition, *independent of moduli*;
   - crack opening δ(x) = 2(1−ν′)σ₀/G·√(L²−x²), hence V = ∫₋Lᴸ δ dx = π(1−ν′)σ₀L²/G.

**Loading decomposition.** Loading 1 (σ_n = −H(t), p = 0 on faces) equals the undrained
elastic crack (instantaneous) plus a diffusion correction of loading-2 type with face
pressure −p_b, where p_b is the undrained pressure induced at the faces. Loading 2
(σ_n = 0, p = H(t) on faces) is that pure diffusion part. The physical pressurized crack
(σ_n = −p_f, p = p_f) is L1 + L2 with magnitude p_f.

## Derivation (brief)

### (a) Loading 1

- **t = 0⁺:** undrained elastic crack: K1_0 = σ₀√(πL), V1_0 = π(1−ν_u)σ₀L²/G.
  The undrained condition ζ = 0 in the fluid-content law gives
  p = −(ν_u−ν)σ_kk/[α(1−2ν)]. For the uniformly loaded Griffith crack (Westergaard
  function Z = σ₀z/√(z²−L²)) the faces have σ_xx = 0, σ_yy = −σ₀, so σ_kk = −σ₀ and
  **p_faces_0 = (ν_u−ν)σ₀/[α(1−2ν)]** > 0 (compression of the skeleton raises p).
  The boundary condition p = 0 on the faces is met through a zero-thickness layer at t = 0⁺.
- **t → ∞:** p → 0 everywhere → drained elastic crack: K1_inf = σ₀√(πL),
  V1_inf = π(1−ν)σ₀L²/G.
- **Short time (ct/L² ≪ 1):** drainage of p_b through the faces is a 1-D boundary-layer
  process of thickness √(ct). Along the faces the field is locally y-dependent only, so
  e_xx = 0 (oedometric). From the constitutive law with e_xx = 0:
  e_yy = (1−2ν)(σ_yy + αp)/[2G(1−ν)], and the mass balance reduces to
  ∂p/∂t = c ∂²p/∂y² with exactly the c quoted in the task (verified in check 4). Hence
  p(y,t) = p_b erf(y/2√ct); the near-face skeleton contracts as pressure drains and pulls
  the faces apart: Δu_y(0⁺,t) = −∫₀^∞ Δe_yy dy = α(1−2ν)p_b√(ct)/[G(1−ν)√π]
  (using ∫₀^∞ erfc(s)ds = 1/√π). With ΔV = 2L·2Δu_y and p_b as above,

  **V ≈ V1_0 + A√(ct/L²),  A = [4(ν_u−ν)]/[(1−ν)√π]**  (units σ₀L²/G; α cancels).

  The crack-tip regions (size ~√(ct)) contribute only O((ct)^{3/4}L^{1/2}) to V, which is
  higher order, so the uniform face contribution above is the leading term.

### (b) Loading 2 (σ_n = 0, p = H(t) on faces)

- **t = 0⁺:** p ≡ 0 in the bulk. Writing total stress σ = σ″ − αpδ, the skeleton sees an
  effective face traction σ″_yy = αp_f, but equilibrium across the (zero-thickness)
  pressure layer, σ″_yy,y = αp_,y, gives σ″_yy just outside the layer = αp_f − αp_f = 0:
  the outer problem feels no loading. **K2_0 = 0** (and V2_0 = 0).
- **t → ∞:** 2-D diffusion exterior to the crack with face value p = p_f maintained has
  no steady state; p(x,t) → p_f pointwise everywhere. Uniform p means uniform effective
  stress σ″ = αp_fδ, zero total stress, uniform strain e_yy = α(1−2ν)p_f/(2G); the
  displacement u_y = e_yy·y vanishes at y = 0±, so the crack faces do not separate:
  **K2_inf = 0, V2_inf = 0**.
  (At intermediate times K2 and V2 are *negative* — the pressurized near-face skeleton
  bulges into the cavity — and relax back to 0. Consistency: for the physical pressurized
  crack L1+L2, K(0⁺) = K(∞) = p_f√(πL) and V grows from π(1−ν_u)p_fL²/G to
  π(1−ν)p_fL²/G, the known undrained→drained behavior.)

### (c) Numerical values, ν = 0.2, ν_u = 0.4, α = 0.89
(K in √(πL), V in L²/G, per unit loading magnitude)

| quantity | value |
|---|---|
| K1_0 | 1 |
| K1_inf | 1 |
| V1_0 = π(1−ν_u) = 0.6π | 1.884955592154 |
| V1_inf = π(1−ν) = 0.8π | 2.513274122872 |
| V1_short_coeff A = 4(0.2)/(0.8√π) = 1/√π | 0.564189583548 |
| p_faces_0 = 0.2/(0.89·0.6) | 0.374531835206 |
| K2_0 | 0 |
| K2_inf | 0 |
| V2_inf | 0 |

## Checks performed

All from `python3 solve.py`; full pasted output in `run_output.txt`. What was compared,
and how each check could have failed:

1. **Crack-volume formula:** numerical quadrature of the COD ellipse
   ∫₋₁¹ 2(1−ν′)√(1−x²)dx vs π(1−ν′), for ν′ = ν and ν_u. A wrong prefactor in V
   (e.g. a missing factor 2, or the wrong modulus dependence) would fail here.
   `[PASS] COD quadrature for V (drained): got 2.51327412193, want 2.51327412287, rel err 3.74e-10`
   `[PASS] COD quadrature for V (undrained): got 1.88495559145, want 1.88495559215, rel err 3.74e-10`
2. **∫₀^∞ erfc(s)ds = 1/√π** by quadrature: validates the boundary-layer integral.
   `[PASS] int_0^inf erfc(s) ds: got 0.56418958355, want 0.564189583548, rel err 4.17e-12`
3. **Oedometric stiffness from the given 2-D constitutive law:** finite-difference slope
   de_yy/dp along the path e_xx = 0, σ_yy fixed, vs α(1−2ν)/[2G(1−ν)]. Would fail if the
   strain–pressure coupling used in the layer analysis were wrong.
   `[PASS] oedometric de_yy/dp: got 0.33374999997, want 0.33375, rel err 8.98e-11`
4. **Diffusivity:** c = κ/(dζ/dp) computed from the constitutive law along the oedometric
   path vs the formula quoted in the task. Would fail if the 1-D reduction of the mass
   balance were inconsistent with the stated c.
   `[PASS] diffusivity c (derived vs given): got 1.87032127464, want 1.87032127444, rel err 1.11e-10`
5. **Independent physics check of the short-time coefficient:** explicit
   finite-difference solution of ∂p/∂t = c∂²p/∂y² on y∈[0,2] (dy = 10⁻³, p(y,0) = p_b,
   p(0,t) = 0), then ΔV(t) = −2Lα(1−2ν)/[G(1−ν)]∫(p−p_b)dy and
   A_num = [ΔV/(L²/G)]/√(ct/L²) vs A = 0.564189583548. Would fail if the sign, the
   factor 4 (= 2 faces × 2L), or the oedometric coupling were wrong.
   `[PASS] ct/L²=0.001: got 0.564253051394 (rel err 1.12e-04); ct/L²=0.003: 0.56421074027
   (3.75e-05); ct/L²=0.01: 0.564195930646 (1.12e-05)` — converging to A as the layer
   thickens relative to dy.
6. **Analytic erf-profile quadrature:** the same ΔV integral with p = p_b erf(y/2√ct)
   vs A. `[PASS] analytic-profile short-time A: got 0.56418958355, want 0.564189583548, rel err 3.47e-12`

`checks finished: 0 failure(s)`.

## Missing or uncertain

- Intermediate-time histories K(t), V(t) were not required and not computed; only the
  endpoints and the leading short-time term are exact. (At short times K1 also picks up a
  near-tip correction of order (ct/L²)^{1/4}, which does not affect any reported item.)
- The t → ∞ limit of loading 2 is pointwise: p = 0 at infinity at every finite time, but
  the diffusion front recedes and p → p_f at every fixed point; the reported long-time
  state is the fully equilibrated uniform-pressure state.
- Sign conventions: tension positive; K > 0 opening mode; V > 0 open; p > 0 fluid
  overpressure. p_faces_0 is positive (pressure rise next to the faces).
- Only the leading term of the short-time expansion is reported; the next term is
  O(ct/L²).

## Questions (with the assumptions used)

1. **Units of `p_faces_0` are not specified** → reported per unit loading magnitude
   (i.e. in units of the applied stress), consistent with K and V being per unit loading.
2. **Loading 2 holds p = p_f on the faces for all t > 0** (maintained reservoir), as
   implied by p = H(t) → in 2-D this eventually uniformizes the pressure, giving the zero
   long-time crack loading reported.
3. **The crack faces are hydraulically connected to the crack cavity**, so p = 0 on the
   faces in loading 1 (as stated) and p = p_f in loading 2.
4. **"Immediately next to the crack faces"** was read as the value of the undrained
   pressure field at y = 0⁺ before diffusion acts; the face value p = 0 itself is
   enforced through a layer of zero thickness at t = 0⁺.
