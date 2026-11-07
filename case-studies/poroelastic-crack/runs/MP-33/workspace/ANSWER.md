# ANSWER — suddenly pressurized crack in a poroelastic solid

**Reproduce command:** `python3 solve.py` (prints every value below, runs 11 self-checks, and
writes `results.json`).

## Model and method, and why

The solid is a linear Biot poroelastic plane-strain medium, quasi-static (no inertia; the stated
constitutive laws have no inertial terms). The two loadings are the two fundamental crack loads —
unit face traction with drained faces (loading 1) and unit face pore pressure with traction-free
faces (loading 2); any real fluid-driven crack is a superposition of them. The response splits
into three regimes, which is why the answers need only elasticity plus a 1-D diffusion layer:

* **t = 0⁺ (undrained):** diffusion has had no time to act, so the bulk responds with ζ = 0.
  Eliminating p from the two given constitutive laws at ζ = 0 gives
  `2G e_ij = σ_ij − ν_u σ_kk δ_ij`, i.e. ordinary elasticity with the **undrained** Poisson
  ratio. (The prescribed face pressure conflicts with the undrained bulk only in a layer of
  zero thickness at t = 0⁺, which cannot affect K or V at t = 0⁺.)
* **t → ∞ (drained, loading 1):** the undrained pressure field is localized and the faces are
  held at p = 0, so p → 0 everywhere; the medium is drained elastic with (G, ν).
* **0 < ct/L² ≪ 1:** diffusion boundary layers of thickness ∼ √(ct) grow from the crack faces;
  this is a local 1-D uniaxial-strain problem.

## Derivation in brief

**(a) Loading 1 (σ_n = −H(t) ≡ unit pressure P = 1 opening the crack, p = 0 on faces).**

*K.* For a pressurized Griffith crack the stress field is the Westergaard solution, independent
of elastic constants (traction BVP). At t = 0⁺ the undrained solid gives K = P√(πL); at t → ∞
the drained solid gives the same. Hence **K1_0 = K1_inf = 1** in units of √(πL) per unit load.

*V.* Plane-strain crack opening displacement under pressure P:
δ(x) = 2(1−ν)P√(L²−x²)/G, so V = ∫δ dx = π(1−ν)PL²/G. With the undrained, resp. drained,
Poisson ratio: **V1_0 = π(1−ν_u)L²/G**, **V1_inf = π(1−ν)L²/G**.

*Pore pressure at the faces at t = 0⁺.* On the faces the elastic solution gives
σ_kk = σ_xx + σ_yy = −2P (Westergaard: σ_xx + σ_yy = 2 Re Z with Re Z = −P on the faces, and
σ_yy = −P). The undrained constraint (ζ = 0) gives p = −(ν_u−ν)σ_kk/[α(1−2ν)], hence
**p_faces_0 = 2(ν_u−ν)/[α(1−2ν)]** (compression ⇒ positive pressure).

*Short-time growth of V.* Fluid is sucked from the solid into the drained crack; a layer of
thickness ∼ √(ct) near each face depressurizes. Locally the layer is in uniaxial strain
(e_xx = 0, uniform along the face since p_faces_0 is constant along it). Then
σ_kk = [σ_yy − α(1−2ν)p]/(1−ν), and mass balance κ∇²p = ∂ζ/∂t reduces to ∂p/∂t = c∂²p/∂y²
with **exactly the stated c** = 2κG(1−ν)(ν_u−ν)/[α²(1−2ν)²(1−ν_u)] — a consistency proof of the
model. With p(y,0) = p0 and p(0,t) = 0: p = p0 erf(y/2√(ct)). The strain is
e_yy = (1−2ν)(σ_yy + αp)/[2G(1−ν)] (Terzaghi effective stress in 1-D), so each face heaves
outward by Δu = (1−2ν)α/[2G(1−ν)] ∫₀^∞(p0−p)dy = (1−2ν)α p0 √(ct)/[G(1−ν)√π]
(using ∫₀^∞ erfc = 1/√π). Both faces over length 2L give ΔV = 4LΔu, i.e.
**V ≈ V1_0 + A√(ct/L²), A = 8(ν_u−ν)/[(1−ν)√π]** (tip zones of size √(ct) contribute only
O(ct/L²), below the order kept).

**(b) Loading 2 (σ_n = 0, p = H(t) on faces).**

*t = 0⁺:* the bulk is undrained with zero applied traction ⇒ no deformation at all; the p = 1
face condition is incompatible with the bulk only in a zero-thickness layer. **K2_0 = 0**
(and V2_0 = 0).

*t → ∞:* a finite crack held at p = 1 in an infinite 2-D medium has no decaying steady state
(exterior Dirichlet problem: the only bounded harmonic function equal to 1 on the crack is
p ≡ 1), so p → 1 pointwise. Then σ → 0 (equilibrium, zero tractions, vanishing at infinity)
and the Biot effective stress σ'_ij = σ_ij + αp δ_ij → αδ_ij, uniform. A uniform effective
stress satisfies the crack-face condition too (σ'_yy = σ_yy + αp = α on the faces), so it
induces **no** crack-tip singularity and no opening: **K2_inf = 0**, **V2_inf = 0** (uniform
free swelling e_ij = α(1−2ν)δ_ij/(2G) leaves the crack invisible). Loading 2 has only a
transient response; note K1_0 + K2_0 = K1_inf + K2_inf = 1, consistent with a real
fluid-pressurized crack (loading 1 + loading 2) having K = √(πL) at both endpoints.

**(c) Numbers for ν = 0.2, ν_u = 0.4, α = 0.89** (K in √(πL), V in L²/G, per unit load):

| quantity | value |
|---|---|
| K1_0 | 1.0 |
| K1_inf | 1.0 |
| V1_0 = π(1−ν_u) | 1.8849555922 |
| V1_inf = π(1−ν) | 2.5132741229 |
| V1_short_coeff = 8(ν_u−ν)/[(1−ν)√π] = 2/√π | 1.1283791671 |
| p_faces_0 = 2(ν_u−ν)/[α(1−2ν)] | 0.7490636704 |
| K2_0 | 0.0 |
| K2_inf | 0.0 |
| V2_inf | 0.0 |

## Checks (run by `python3 solve.py`; outputs pasted verbatim)

1. **Undrained reduction** — compared the coefficient of σ_kk in 2G e_ij after eliminating p at
   ζ=0 against −ν_u: `[check] undrained reduction coefficient of sigma_kk vs -nu_u: got -0.4,
   expected -0.4 -> PASS`. Fails if the two given constitutive laws were mutually inconsistent.
2. **Diffusivity** — compared the c implied by my 1-D uniaxial-strain mass balance against the
   stated formula: `[check] 1-D boundary-layer diffusivity vs stated c: got 0.000579799595075,
   expected 0.000579799595075 -> PASS`. Fails if the layer mechanical constraint (e_xx = 0)
   were the wrong one.
3. **Crack-volume formula** — numerical quadrature of the COD vs π(1−ν)PL²/G:
   `[check] crack volume quadrature vs pi(1-nu)PL^2/G: got 2.51327412287, expected
   2.51327412287 -> PASS`. Fails if the COD expression or its integral (πL²/2) were wrong.
4. **∫erfc** — `[check] integral_0^inf erfc(s) ds vs 1/sqrt(pi): got 0.564189583548, expected
   0.564189583548 -> PASS`. Fails if the heave integral were wrong.
5. **Westergaard face stress** — max |σ_kk + 2P| over 21 face points:
   `[check] max |sigma_kk + 2P| on faces (Westergaard): got 2.23774607733e-08, expected 0 ->
   PASS`, and p_faces_0 recomputed from that trace: `[check] p_faces_0 from Westergaard trace:
   got 0.749063670412, expected 0.749063670412 -> PASS`. Fails if σ_kk on the faces differed
   from −2P, which would change p_faces_0 and A.
6. **Independent finite-difference diffusion simulation** (explicit FTCS, p(y,0)=p0, p(0,t)=0,
   five grids Ny = 1000…16000) of the face layer; dimensionless slope of ΔV vs √(ct/L²) at
   ct/L² = 5×10⁻⁴: `1.106272, 1.117262, 1.122805, 1.125588, 1.126983`, successive-difference
   ratios `0.504, 0.502, 0.501` (first order in dy because of the discontinuous initial
   condition), Richardson-extrapolated `[check] Richardson-extrapolated slope: 1.128377
   (analytic A = 1.128379)` -> PASS. Fails if any element of the layer physics (diffusivity,
   heave law, both-faces factor 4L, erfc profile) were wrong.
7. **erf similarity profile** — quadrature of ∫p0 erfc(y/2√(ct))dy vs 2p0√(ct/π) at
   ct = 1e-4, 1e-3, 1e-2, e.g. `got 0.0845227840521, expected 0.0845227840521 -> PASS`.
   Fails if the similarity solution were wrong.

Final line of the run: `11 checks, 0 failed`.

## What is missing / uncertain

* Intermediate-time behavior (0 < t < ∞) was not required and not solved; K1(t) may deviate
  transiently between its equal endpoints, and V1(t) is monotone between V1_0 and V1_inf only
  by physical expectation, not proven here.
* The short-time law is leading order only; O(ct/L²) corrections from the crack-tip boundary
  layers were not computed.
* Loading 2's transient (which is nonzero at intermediate times) was not computed; only the
  endpoints were required.

## Questions (with assumptions used)

* **Loading magnitude/units:** taken as unit (per unit loading magnitude), as stated for the
  report; all results are linear in the load.
* **Sign convention:** σ_n = −H(t) is taken as unit pressure *opening* the crack (K positive).
* **Loading 2 at t → ∞:** assumed the p = 1 condition is maintained forever; in an infinite 2-D
  medium p then tends to 1 at every finite point (no bounded steady state exists), giving the
  stated zero long-time crack response. If instead a far-field pressure sink at finite distance
  were assumed, K2_inf and V2_inf would not be zero — but that contradicts "vanish at
  infinity".
* **Quasi-static response:** "just after loading" is read as the undrained elastic limit,
  ignoring inertial waves (standard for Biot quasi-static problems).
* **Fully permeable crack faces** assumed, since p is prescribed on them.
