# Suddenly pressurized crack in a Biot poroelastic solid: solution

All derivations below were done from the problem statement only. Numerical work (Python,
`code/`) is used where no closed form was found, and every number is cross-checked
against at least one independent identity or limit (Section 7).

Units used in the summary: K in units of (loading)·√(πL), V in units of (loading)·L²/G,
time as the dimensionless group τ = ct/L².

## 0. Results

Constants used throughout (plane strain, in-plane traces σ_kk = σ_xx + σ_yy):

| symbol | definition | value (ν = 0.2, ν_u = 0.4, α = 0.89) |
|---|---|---|
| η | α(1−2ν)/[2(1−ν)] (poroelastic stress coefficient) | 0.33375 |
| β | (ν_u−ν)/[α(1−2ν)] (undrained: p = −β σ_kk) | 0.374532 |
| γ | (ν_u−ν)/(1−ν) = 2βη | 0.25 |
| S | α²(1−2ν)²(1−ν_u)/[2G(ν_u−ν)(1−ν)] (1-D storage; c = κ/S) | 0.534668/G |

| quantity | general result | value |
|---|---|---|
| K1(0⁺) | a(γ)·σ√(πL); a defined by the drained-tip inner problem (Sec. 5); a = 1 − 3γ/4 + O(γ²); exact a(γ) computed numerically | **0.80759** |
| K1(∞) | σ√(πL) | **1** |
| V1(0⁺) | π(1−ν_u)σL²/G | **1.884956** |
| V1(∞) | π(1−ν)σL²/G | **2.513274** |
| A (V1 ≈ V1(0⁺) + A√τ) | A = (16/√π)γ(1−ν_u) [face layers, exact] + A_tip [crack-tip regions]; A_tip = 4√π γ(1−ν_u)[1+O(γ)] | **2.440** (1.3541 + 1.0862) |
| p next to faces at 0⁺ (loading 1) | 2(ν_u−ν)σ/[α(1−2ν)] = 2βσ | **0.749064** |
| K2(0⁺) | 0 (K2 grows as −C τ^{1/4}; C = 0.2599 here) | **0** |
| K2(∞) | −η p√(πL) = −α(1−2ν)p√(πL)/[2(1−ν)] | **−0.33375** |
| V2(∞) | −πη(1−ν)pL²/G = −πα(1−2ν)pL²/(2G) | **−0.838805** |

Two points need care when these numbers are compared with a simpler analysis:

1. **K1(0⁺) is not the undrained value.** At t = 0 exactly the response is undrained and
   K = σ√(πL). For any t > 0 a region of size ~√(ct) around each tip is drained, and the tip
   K is a(γ)σ√(πL) = 0.8076σ√(πL) no matter how small t is (the tip problem is
   self-similar). So K jumps at t = 0 and then rises:
   K1 ≈ [0.8076 + 0.195 τ^{1/4} + …]σ√(πL). The limit t → 0⁺ asked for is 0.8076.
2. **The √t growth of V1 has two sources of the same order.** The drained face boundary
   layers give A_face = (16/√π)γ(1−ν_u) = 1.354. A naive "uniaxial swelling" estimate,
   which ignores how the outer body responds to the layer, gives (8/√π)γ = 1.128 and is
   incorrect. The crack-tip regions add A_tip = 1.086. So A = 2.440. Using only the
   first-order (weak-coupling) results gives the closed form
   A ≈ 4(4+π)γ(1−ν_u)/√π = 2.418, which is 0.9 % low.

The long-time limits are approached slowly, and in loading 2 only logarithmically (2-D diffusion):
K2(t) ≈ −η[1 − 1/(ln(16τ) − γ_E)], V2(t) ≈ −πη(1−ν)[1 − 1/(ln(16τ) − γ_E)]
(γ_E = Euler's constant).

---

## 1. Reformulations of the constitutive law

The given law is Biot's law already reduced to plane strain (e_zz = 0): σ_kk = σ_xx + σ_yy
(checked by eliminating σ_zz from the 3-D Rice–Cleary law). Useful consequences:

* **Undrained (ζ = 0):** p = −β σ_kk, and 2G e_ij = σ_ij − ν_u σ_kk δ_ij. This is plane-strain elasticity with ν_u.
* **Eigenstrain form:** eliminating p gives

  2G e_ij = σ_ij − ν_u σ_kk δ_ij + 2G β ζ δ_ij.

  So the poroelastic body is an undrained elastic body carrying the isotropic in-plane
  eigenstrain βζ.
* **Diffusion of ζ:** equilibrium and compatibility give ∇²(σ_kk + 2ηp) = 0. Hence
  ζ = S p + α(1−2ν)(σ_kk + 2ηp)/(2G) with the second term harmonic, and
  ∂ζ/∂t = κ∇²p becomes **∂ζ/∂t = c∇²ζ** exactly, with c = κ/S equal to the given diffusivity.
* **1-D consolidation layer** (σ_yy fixed, e_xx fixed): Δσ_xx = −2ηΔp, Δe_yy = ηΔp/G,
  Δζ = SΔp, and Δp obeys ∂_tΔp = c∂_yyΔp.

**Reciprocal identity for the volume.** Apply Betti's theorem in the undrained elastic body
between the actual state (applied face traction plus eigenstrain βζ) and the auxiliary
state "unit pressure on the faces, undrained". This gives

  V(t) = V_el + ∫ σ^a_kk β ζ dA = V_el − ∫ p^a ζ dA,                                  (R)

where V_el is the undrained elastic volume for the applied face traction (zero in loading 2),
σ^a is the undrained stress field of the unit-pressure Griffith crack, and
p^a = −βσ^a_kk = −2β Re[z/√(z²−L²) − 1] is its undrained pore pressure.
p^a equals 2β on the faces and is −2βK cos(θ/2)/√(2πr) + 2β near a tip, with K = √(πL).
The identity is exact for all t.

## 2. Undrained (t = 0⁺) and drained (t → ∞) Griffith solutions

A uniformly pressurized Griffith crack has σ_xx = σ_yy = −σ on the faces, K = σ√(πL) for any
Poisson ratio, and opening w = 2(1−ν')σ√(L²−x²)/G, so V = π(1−ν')σL²/G.

* **Loading 1, t = 0⁺.** Away from vanishing boundary layers the response is undrained (ν' = ν_u):
  **V1(0⁺) = π(1−ν_u)σL²/G**. The induced pore pressure next to the faces follows from
  σ_kk = −2σ: **p = 2βσ = 2(ν_u−ν)σ/[α(1−2ν)]** (equivalently 2B(1+ν_u)σ/3 with
  Skempton's B). This is incompatible with p = 0 on the faces, so it starts to diffuse
  into the crack.
* **Loading 1, t → ∞.** The pore pressure decays to 0. The transient has finite total fluid
  content, so in 2-D its effect decays like 1/t. The state is the drained Griffith solution:
  **K1(∞) = σ√(πL), V1(∞) = π(1−ν)σL²/G**.

## 3. Loading 2 at long times: the poroelastic back-stress

The crack keeps p = 1 on its faces. In 2-D the pressurized zone keeps growing (radius
~√(ct)), and near the crack p → 1 (the deficit decays like 1/ln t). Use the thermoelastic
analogy α(1−2ν)p ↔ Eα_T T. A large disc of radius R with uniform p in an infinite plane,
stress-free at infinity, has uniform interior stress σ_xx = σ_yy = −ηp, independent of R.
Uniform p by itself produces no stress and no opening (free swelling). So near the crack the
stress is an isotropic compression −η. The crack faces must stay free of total normal
traction, so the correction problem is a Griffith crack whose faces carry the tension +η.
Hence

  **K2(∞) = −η p√(πL),  V2(∞) = −π(1−ν)η pL²/G = −πα(1−2ν)pL²/(2G).**

Check with (R): at long times near the crack, ζ2 = S − α(1−2ν)ησ^a_kk/(2G). The Griffith
integrals ∫σ^a_kk dA = −πL² and ∫(σ^a_kk)² dA = 2πL² were verified numerically
(`code/griffith_integrals.py`). Then V2 = −πL²[η(1−ν_u) + η(ν_u−ν)]/G = −πη(1−ν)L²/G. ✓

For the combined fluid-filled crack (σ_n = −p_f, p = p_f), this is the familiar back-stress:
K → (1−η)p_f√(πL).

**Rate of approach (numerical, Sec. 6).** The Laplace transforms satisfy, to 6 digits,
sK̃2 = −η[1 − 1/(ln(1/s) + 4 ln 2 − 2γ_E)], and the same form holds for V2. Here
4 ln 2 − 2γ_E is the logarithmic-capacity constant of a segment of length 2L. In the time
domain this gives K2 ≈ −η[1 − 1/(ln(16τ) − γ_E)].

## 4. Small times, loading 1: growth of the crack volume

For ct ≪ L² the field has three parts:
* the undrained outer solution of Sec. 2;
* drained boundary layers of thickness ~√(ct) along the faces;
* tip regions of size ~√(ct).

Outside these regions ζ = 0 exactly: in the Laplace domain the diffusive part decays like
exp(−d√(s/c)) with distance d. So (R) needs ζ only inside the layers and tip regions.

**4.1 Face layers (exact to O(√t)).** In the layer, Δp = p − 2βσ = −2βσ erfc(y/2√(ct)),
Δζ = SΔp, and p^a = 2β on the faces. Over 2 faces × length 2L,

  ΔV_face = 2β·S·2β·(2√(ct/π))·4L = 32β²S L√(ct/π) = (16/√π)·γ(1−ν_u)·(σL²/G)·√(ct/L²),

using β²SG = γ(1−ν_u)/2. A direct "uniaxial swelling" estimate (wall shrinkage ∫Δe_yy dy only)
gives (8/√π)γ instead. The difference, (8/√π)γ(1−2ν_u), is the elastic response of the
outer body to the tensile membrane force N = −2η∫Δp dy carried by the constrained,
contracting layer; that force is transmitted to the body at the tips. (R) contains this
response automatically. Numerically, loading 2 (Sec. 6) reproduces the coefficient
−(16/√π)η(1−ν_u) of its own face-layer law V2 ≈ −(16/√π)η(1−ν_u)√τ to about 0.1 %, which is
the extrapolation accuracy.

**4.2 Tip regions (the same order).** Near a tip the undrained pore pressure is
p_K = −2βK cos(θ/2)/√(2πr) (suction). It is harmonic, vanishes on the faces, and has ζ = 0,
so it formally satisfies all field equations. But it needs an infinite fluid supply at the tip:
the face flux ~ r^{-3/2} is non-integrable and the dissipation ∫|∇p|² diverges. For t > 0 it
is therefore regularized within r ~ √(ct). At first order in the coupling (mechanical field
frozen), p obeys ∂_t p = c∇²p with p(0⁺) = p_K, p = 0 on the faces, and regularity at the tip.
Expanding in the angular mode cos(θ/2) (a Hankel transform of order ½) gives the exact solution

  p = p_K erf(r/2√(ct)),   ζ = S(p − p_K) = −S p_K erfc(r/2√(ct)).

Fluid is drawn from the crack into the suction zone around the tip, and the swelling there
opens the crack. From (R), for two tips with K² = πσ²L:

  ΔV_tip = 2S∫p_K² erfc(r/2√(ct)) dA = 8β²SK²√(ct/π) = 4√π·γ(1−ν_u)·(σL²/G)·√(ct/L²)  [first order].

With full coupling, the tip coefficient is A_tip = (2/√π)(1−ν_u)ĉ(γ). The numerical
solution (Sec. 6) gives ĉ(γ) = 2πγ + 0.434γ² + 0.37γ³ (fit, |error| < 5·10⁻⁴ for γ ≤ 0.5)
and ĉ(0.25) = 1.6043, so A_tip = 1.0862.

**4.3 Result.**

  V1 ≈ π(1−ν_u)σL²/G + A·(σL²/G)·√(ct/L²),
  A = (16/√π)γ(1−ν_u) + A_tip ≈ (4/√π)(4+π)·(ν_u−ν)(1−ν_u)/(1−ν)  (first order in γ),

which evaluates to A = 1.3541 + 1.0862 = **2.440** (first-order closed form 2.4175). The next
term is O((ct/L²)^{3/4}). A direct extrapolation of the full Laplace solution of loading 1 gives
2.4400, which is consistent within the extrapolation error of about 0.0003.

For loading 2 the leading term is pure face layer: V2 ≈ −(16/√π)η(1−ν_u)(pL²/G)√τ
= −1.8077√τ. Its tip contributions start at O(τ^{3/4}).

## 5. Stress intensity factor at t = 0⁺

**5.1 Inner problem.** For ct ≪ L², each tip sees a semi-infinite crack:
* far field (r ≫ √(ct)): the undrained K-field with K∞ = σ√(πL), plus O(1) uniform fields;
* faces: traction-free and p = 0;
* near field (r ≪ √(ct)): drained.

The only length is √(ct), so K_tip = a·K∞ is constant in time. The O(1) uniform part (the face
mismatch 2β) is −2β times the loading-2 inner problem. It adds only O((ct)^{1/4}).

In scaled variables the Laplace-domain kernel of this problem depends on a single parameter:

  [σ̂', p̂'] = [[−|k|(1−γ·2k²/(m(m+|k|))), −√(2γ)·|k|/(m(m+|k|))], [+√(2γ)·|k|/(m(m+|k|)), 1/m]]·[Û', Q̂'].

So **a = a(γ)**. This was confirmed numerically: four (ν, ν_u, α) sets with γ = 0.25 give
identical results to 8 digits.

**5.2 Laplace-domain J-integral (includes a fluid term).** The transformed equations are the
Euler–Lagrange equations of 𝓛 = W_d(ε) − αpe − p²/(2M) − (κ/2s)|∇p|², where
1/M = α²(1−2ν)(1−2ν_u)/[2G(ν_u−ν)]. So
P_xj = 𝓛δ_xj − σ_ij u_i,x + (κ/s) p,_j p,_x is divergence-free. Its crack-face terms vanish
because the faces are traction-free and p is constant along them. The tip field is the drained
K-field plus p ≈ A_p√r cos(θ/2). Therefore

  (1−ν)K̃_tip²/(2G) + πκÃ_p²/(4s) = (1−ν_u)K̃∞²/(2G),

and a² = (1−ν_u)/(1−ν) − (fluid term) < (1−ν_u)/(1−ν). The value √((1−ν_u)/(1−ν)) = 0.866 is
only an upper bound.

**5.3 Weak coupling.** From the erf solution, A_p = −2βK/(π√(2ct)) and Ã_p = −2βK/√(2πcs). The fluid
term is then exactly γ/2 of the far-field J, so a² = (1−γ)(1−γ/2) + O(γ²):

  **a = 1 − (3/4)γ + O(γ²)**.

**5.4 Exact value.** The full inner problem is a 2×2 matrix Wiener–Hopf problem, and no closed
factorization was found. a(γ) was therefore obtained as lim_{s→∞} of the finite-crack Laplace
solution for the loading σ_n = −1, p = 2β. That loading has no face layer, so the limit converges
like s^{-1/2}, which makes the extrapolation clean.

| γ | 0.01 | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.40 | 0.50 |
|---|---|---|---|---|---|---|---|---|---|
| a(γ) | 0.99249 | 0.96232 | 0.92428 | 0.88583 | 0.84695 | **0.80759** | 0.76772 | 0.68621 | 0.60185 |
| ĉ(γ) | 0.0629 | 0.3154 | 0.6332 | 0.9538 | 1.2773 | **1.6043** | 1.9341 | 2.6062 | 3.2971 |

Fit: a ≈ 1 − 0.75γ − 0.0634γ² − 0.0582γ³ (|error| < 10⁻⁴ for 0 < γ ≤ 0.5). The values satisfy
(1−γ)^{3/4} < a < √((1−γ)(1−γ/2)); for γ ≤ 0.25 both of these simple forms are within about 0.3 % of a.

For the given data (γ = 0.25): **K1(0⁺) = 0.80759 σ√(πL)**. The approach to this value is slow:
K1(t) ≈ [0.8076 + 0.1947τ^{1/4} + …]σ√(πL); for example K1 = 0.827 at τ = 10⁻⁴ (Talbot inversion).

**5.5 Loading 2: K2(0⁺) = 0.** The inner problem is driven only by the face pressure, whose
layer has thickness √(ct). The far-field J of that problem is 2S∫p²dy = O(√(ct)), so
K2 = O((ct)^{1/4}) → 0. Numerically K2 ≈ −0.2599·τ^{1/4}·p√(πL), and its Laplace transform
satisfies s^{1/4}·sK̃2 → −0.23554. K2 is negative because the swelling walls close the crack.

## 6. Numerical solution of the full problem (Laplace domain)

**Formulation.** Laplace transform in t (parameter s); symmetric about y = 0, so solve in the
half-plane y > 0. The general solution there is
* an undrained elastic field (Papkovich–Neuber, ν_u), with ζ = 0 and p harmonic, which solves
  the transformed equations exactly; plus
* one irrotational diffusive mode: u = ∇φ, ∇²φ = (s/c)φ, p = (G/η)(s/c)φ,
  σ = 2G(∇∇φ − ∇²φ I).

On y = 0, impose σ_xy = 0 and write U = u_y(x, 0⁺) (half the opening) and Q = q_y(x, 0⁺)
(flux from the crack into the solid). The Fourier symbols, derived with sympy in
`code/symbols.py` (fields ∝ e^{ikx}, m = √(k² + s/c)), are

  σ̂_yy = A11Û + A12Q̂,   p̂ = A21Û + A22Q̂,
  A11 = −(G|k|/(1−ν_u))·[1 − 2γk²/(m(m+|k|))],   A22 = 1/(κm),
  A21 = (2Gβ/(1−ν_u))·|k|(m−|k|)/m,   A12 = −A21/s.

They reduce to undrained elasticity as s → ∞ and to drained elasticity as s → 0 or |k| → ∞.

Mixed conditions:
* |x| < L: σ_yy and p prescribed (loading 1: −1/s, 0; loading 2: 0, 1/s);
* |x| > L: U = Q = 0.

**Galerkin method** (L = G = c = 1):
* Basis: U = Σa_n√(1−x²)U_n(x) and Q = Σb_nT_n(x)/√(1−x²), with even n. These carry the
  r^{1/2} opening and r^{-1/2} tip flux.
* Projection: the two equations are projected on the same functions; the transforms are
  π(n+1)(−i)ⁿJ_{n+1}(k)/k and π(−i)ⁿJ_n(k).
* Matrix entries are ∫_0^∞(symbol)·J_μJ_ν·k^{-j}dk:
  - the drained/large-k parts are done in closed form (Weber–Schafheitlin);
  - the remainders by Gauss–Legendre on [0, k_max], k_max = 4(2N+1) + 40√|s| + 200, with
    graded panels near k = 0 when |s| < 1;
  - the tail k > k_max with the non-oscillatory Debye form of J_μJ_ν;
  - the (0,0) flux entry exactly, (2/πκ)∫_0^{π/2} I₀K₀(√s cos θ)dθ, via Neumann's product formula.
* Outputs: K/√(πL) = GΣ(n+1)a_n/(1−ν) (the tip is drained) and V = πa_0.
* Time histories: fixed-Talbot inversion (M = 12, complex s).

Code: `code/crack_solver2.py` (solver), `code/talbot.py`, `code/time_domain.py`,
`code/final_extrapolation.py`, `code/a_of_gamma.py`, `code/fit_gamma.py`,
`code/small_s.py`, `code/energy_check.py`, `code/check_gamma_scaling.py`,
`code/griffith_integrals.py`, `code/check_erf_solution.py`, `code/run_large_s.py`,
`code/write_results.py`. Earlier diagnostics are kept as `code/crack_solver.py` (first
version, insufficient k-tail) and `code/conv_test*.py`, `code/diag*.py`.

**Time histories** (ν = 0.2, ν_u = 0.4, α = 0.89; units as in Sec. 0):

| τ = ct/L² | K1 | V1 | K2 | V2 |
|---|---|---|---|---|
| 0⁺ (limit) | 0.8076 | 1.8850 | 0 | 0 |
| 0.0001 | 0.8273 | 1.9065 | -0.0260 | -0.0162 |
| 0.001 | 0.8430 | 1.9469 | -0.0460 | -0.0473 |
| 0.003 | 0.8544 | 1.9849 | -0.0604 | -0.0772 |
| 0.01 | 0.8711 | 2.0490 | -0.0811 | -0.1293 |
| 0.03 | 0.8912 | 2.1336 | -0.1057 | -0.2007 |
| 0.1 | 0.9197 | 2.2535 | -0.1399 | -0.3093 |
| 0.3 | 0.9505 | 2.3682 | -0.1774 | -0.4269 |
| 1 | 0.9792 | 2.4574 | -0.2181 | -0.5435 |
| 3 | 0.9925 | 2.4940 | -0.2459 | -0.6169 |
| 10 | 0.9978 | 2.5076 | -0.2656 | -0.6673 |
| 100 | 0.9998 | 2.5127 | -0.2866 | -0.7203 |
| 1000 | 1.0000 | 2.5132 | -0.2979 | -0.7487 |
| 10000 | 1.0000 | 2.5133 | -0.3049 | -0.7663 |
| ∞ (limit) | 1 | 2.5133 | -0.3338 | -0.8388 |

All four histories are monotone. K1 and V1 are close to drained by τ ≈ 10. Loading 2 creeps
toward its limits like 1/ln τ.

## 7. Checks performed

1. **Limits.** As s → 0: sK̃1 → 1 and sṼ1 → π(1−ν) to 10⁻⁹. As s → ∞: sṼ1 → π(1−ν_u).
2. **Energy identity (Laplace domain).** J computed at the tip from the Galerkin coefficients,
   J = (1−ν)πK̃²/2 + π(Σb_n)²/(2κs), equals (1/4s)∂Ṽ/∂L = Ṽ/s + ½dṼ/ds.
   * The agreement is 10⁻⁶ relative at s = 0.1, 1, 10, 100.
   * The fluid tip term is essential: it is 4.6 % of J at s = 100.
   * This independently confirms both the solver and the fluid term in the J-integral of Sec. 5.2.
3. **Convergence.** Results are unchanged to ≤ 10⁻⁷ (in sK̃) when the number of modes is
   increased (N = 20…190), k_max is doubled, Gauss points are doubled, or θ points are increased.
   * A first version without the Debye tail showed 10⁻⁴ noise in K; it is kept only as a record.
   * Small s needed graded k-panels. An ungraded version gave a wrong long-time K2, and that bug
     was found because the result contradicted Sec. 3.
4. **Single-parameter reduction.** Four parameter sets with γ = 0.25 give identical
   W = s(K̃1 + 2βK̃2) and V_W/(1−ν_u) at s = 10² and 10⁴.
5. **First-order theory against numerics** (small γ):
   * (1−a)/γ = 0.75069 at γ = 0.01 and → 0.7500 as γ → 0 (theory 3/4);
   * ĉ/(2πγ) = 1.0008 at γ = 0.01 (theory 1).
6. **Face-layer law.** The extrapolated lim √s·sṼ2 = −1.6019 (data up to s = 10⁷) matches the exact
   −8η(1−ν_u) = −1.6020.
7. **Long-time fits.** Fitting X∞ + C/(ln(1/s) + D) to s = 10⁻⁸…10⁻¹² gives:
   * X∞ = −0.333750 for K2 and −0.838805 for V2;
   * D = 1.6181 = 4 ln 2 − 2γ_E, C = −X∞ in both.
8. **Griffith integrals.** ∫σ_kk dA = −πL² and ∫σ_kk² dA = 2πL², used in the reciprocal check of V2(∞).
9. **First-order tip solution (symbolic).** sympy confirms that
   p = r^{-1/2}cos(θ/2)·erf(r/2√(ct)) solves ∂_t p = c∇²p, vanishes on θ = ±π, and behaves like
   √r/√(πct) at the tip (`code/check_erf_solution.py`).
10. **Time domain.** The Talbot histories agree with the small-time laws
   (K1 = 0.8271 predicted vs 0.8273 at τ = 10⁻⁴; K2 = −0.0260 vs −0.0260) and with the long-time
   law (K2 = −0.3045 predicted vs −0.3049 at τ = 10⁴).

## 8. Uncertainties and caveats

* **No closed form for a(γ) or A_tip.** a(γ) (K1 at 0⁺) and the tip part of A were not obtained
  in closed form. The inner problem is a genuinely coupled 2×2 Wiener–Hopf problem.
  * Exact here: their first-order (weak-coupling) expressions and the J-integral bound.
  * Numerical: their values at γ = 0.25.
  * Estimated uncertainties: a = 0.807593 ± 0.000005; A_tip = 1.0862 ± 0.0003, so A = 2.4401 ± 0.0003.
    The slow s^{-1/4} corrections in V limit the precision of A.
* **Meaning of "t = 0⁺" for K.** If it is read as the instantaneous undrained response (t = 0
  exactly), the answer is K1 = σ√(πL), the outer undrained value. The limit lim_{t↓0}K1(t) is
  0.8076σ√(πL). That limit is reported because the near-tip region is drained for every t > 0.
  V is continuous, so V1(0⁺) is the same under either reading.
* **No contact constraint.** In loading 2 the linear model lets the crack "close" (V2 < 0);
  this is meaningful only in superposition with loading 1 (fluid-filled crack).
* **Long-time limits.** K2(∞) and V2(∞) are reached only logarithmically in time.
