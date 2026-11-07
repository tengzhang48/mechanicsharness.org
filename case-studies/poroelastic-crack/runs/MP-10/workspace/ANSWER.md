# ANSWER — suddenly loaded / pressurized crack in a poroelastic solid

**Model & method.** Linear Biot poroelasticity, plane strain, as given. The two limits in time are
exact elastic analogies, and the loading-2 steady state admits a closed-form exact solution:

* **t = 0⁺ (undrained):** every field with `ζ = 0`; the response is the *elastic* Griffith-crack
  problem with Poisson's ratio `ν_u` (same G). For a Griffith crack of half-length L under uniform
  face loading the stress field is independent of the elastic constants:
  `K_I = σ√(πL)`, `V = π(1−ν)σL²/G`, and on the crack faces `σ_xx = σ_yy = σ_n`, `σ_kk = 2σ_n`
  (from the Westergaard function `Z = σ_n(1 − z/√(z²−L²))`).
* **t → ∞ (drained/steady):** p solves Laplace's equation with the face boundary values; with
  `p = 0` on the faces (loading 1) p ≡ 0 and the problem is elastic with ν.
* **Loading 2 steady state:** `p_ss = (2/π)Re[arcsin(L/z)]` (harmonic, =1 on the slit, →0 at ∞).
  I verify that **σ_ij ≡ 0** with strain `e_ij = α(1−2ν)p_ss/(2G) δ_ij` satisfies every equation
  and boundary condition (the strain is compatible because p_ss is harmonic — the poroelastic
  analogue of the stress-free steady thermal field). Hence `K2_inf = 0` and the crack volume is the
  pure swelling displacement `u = c(Re f, Im f)`, `c = α(1−2ν)/(2G)`, `f′ = (2/π)arcsin(L/z)`,
  with tip closure `u_y(±L)=0` giving `V2_inf = 4cL²/π = 2α(1−2ν)L²/(πG)`.
* **Short-time V1:** a 1-D consolidation skin of thickness ~√(ct) forms on the faces (the crack is
  otherwise undrained). On the faces the undrained state is `σ_yy = σ_xx = −1`, `p* = 2(ν_u−ν)/(α(1−2ν))`
  (Skempton: `p_u = −(ν_u−ν)σ_kk/(α(1−2ν))`, `σ_kk = −2`). With `p(0)=0` the skin obeys
  `p = p* erf(y/2√(ct))` with exactly the diffusivity c given (checked), and the face displacement
  grows by `Δu = α(1−2ν)p*√(ct)/(G(1−ν)√π)` per face (coefficient independent of the lateral
  strain constraint). Tip zones contribute only O((ct/L²)^{3/4}). Hence
  **`V1 ≈ V1_0 + [8(ν_u−ν)/((1−ν)√π)]·(L²/G)·√(ct/L²)`** (α cancels).

## Answers (c): ν = 0.2, ν_u = 0.4, α = 0.89  (K in units √(πL), V in units L²/G, per unit |load|)

| quantity | value | formula (general) |
|---|---|---|
| K1_0  | **1.000000** | σ√(πL) — undrained limit; K of a Griffith crack is independent of elastic constants |
| K1_inf| **1.000000** | σ√(πL) — drained limit (same reason) |
| V1_0  | **1.884956** | π(1−ν_u)σL²/G |
| V1_inf| **2.513274** | π(1−ν)σL²/G |
| V1_short_coeff A | **1.128379** | 8(ν_u−ν)/((1−ν)√π) |
| p_faces_0 | **0.749064** | 2(ν_u−ν)/(α(1−2ν))·σ  (= 2B/3, B = 1.1236) |
| K2_0  | **0.000000** | no mechanical load at 0⁺; pore-pressure layer has zero thickness |
| K2_inf| **0.000000** | steady state has total stress ≡ 0 (exact solution above) |
| V2_inf| **0.339955** | 2α(1−2ν)L²/(πG) |

Note (not requested): at intermediate times things are *not* zero — K2(t) rises from 0 and returns
to 0, and V2(t) first dips (short-time coefficient −4α(1−2ν)/((1−ν)√π) = −1.5064) before
approaching +0.3400. The crack volume in loading 1 grows monotonically 1.885 → 2.513.

## Reproduce

```
python3 solve.py
```
prints all values above and the checks below, and writes `results.json`.

## Checks performed (all pass; outputs pasted from `python3 solve.py`)

1. **Elastic crack K and face stresses** (Westergaard `Z=−(1−z/√(z²−L²))`), could fail if the
   stress function or K-extraction were wrong:
   `x-L=1e-06: K_num/sqrt(pi L) = 0.998587, bias-corrected = 1.00000075 (-> 1.0)`;
   `sigma_xx on faces = ['-1.000000','-1.000000','-1.000000','-1.000000']` → σ_kk(faces) = −2 → p_faces_0 = 0.74906.
2. **Elastic crack volume from COD** `δ(x)=2(1−ν)√(L²−x²)/G` (could fail if COD wrong):
   `V1_0 (nu_u): numeric = 1.8849555922, formula pi(1-nu) = 1.8849555922`;
   `V1_inf (nu): numeric = 2.5132741229, formula = 2.5132741229`.
3. **p_ss = (2/π)Re[arcsin(L/z)]** (could fail if not harmonic / wrong BCs, which would kill the
   σ≡0 solution): `Laplacian(p_ss) at 4 points = ['-4.44e-08','-2.22e-08','2.22e-08','0.00e+00']`;
   `p_ss on slit x=0.3,0.6,0.9: [1.0,1.0,1.0]`; decays to 6.3e-2 at (5,5), 2.6e-2 at (20,10).
4. **V2_inf** (could fail if the swelling-displacement construction were wrong):
   `I_num = 0.5000000000 vs L^2/2 = 0.5` (∫₀^L x·arccosh(L/x)dx);
   `V2_inf numeric (2*int u_y dx) = 0.33995529 vs formula 0.33995496`.
5. **Skin integral for A** (could fail if the 1-D consolidation bookkeeping were wrong):
   `int erfc(y/(2 sqrt(ct))) dy: numeric = 0.01128379, 2*sqrt(ct)/sqrt(pi) = 0.01128379`
   → A = 4α(1−2ν)p*/((1−ν)√π) = 1.12837917.
6. **Skempton consistency:** B = 3(ν_u−ν)/(α(1−2ν)) = 1.12359551, p_u(face) = 2B/3 = 0.74906.
7. **Independent symbolic re-derivation of the whole A-chain from the problem's constitutive
   laws** (`python3 check2.py`, sympy; could fail if any hand algebra slip in p*, the skin
   strain coefficient, the diffusivity, or the volume bookkeeping):
   ```
   (1) p* (undrained face pressure) = 2*sigma*(nu - nu_u)/(alpha*(2*nu - 1))
       check p* = 2B*sigma/3 : True
   (2) d e_yy/d p (skin, any fixed e_xx) = alpha*(2*nu - 1)/(2*G*(nu - 1))
       == alpha(1-2nu)/(2G(1-nu)): True
   (3) skin diffusivity == given c: True
   (4) du_face/(sqrt(c t)) = 2*sigma*sqrt(t)*(nu - nu_u)/(sqrt(pi)*G*(nu - 1))
   (5) A = 8*sigma*(nu - nu_u)/(sqrt(pi)*(nu - 1))  == 8 sigma (nu_u-nu)/((1-nu)sqrt(pi)): True
       numeric A = 1.12837916710
       numeric p* = 0.749063670412
       V2_inf = 2*alpha*(1-2nu)/pi = 0.339954958444
   ```

## Questions / assumptions

* "per unit loading magnitude" read as: values per unit |σ_n| = 1 (loading 1) and p = 1 (loading 2),
  with the sign convention that compressive σ_n (internal-pressure-type loading) opens the crack
  (V > 0), consistent with "positive for an open crack".
* K is the mode-I intensity of the **total** stress, σ_yy → K/√(2πr) ahead of the tip.
* Both crack faces carry identical conditions (symmetric loading), faces fully connected to the
  crack pore pressure (p = 0 / p = 1 as stated).
* K2_0 = 0 exactly (no instantaneous mechanical response to a pore-pressure step); the problem's
  request of only the two time-endpoints of K2 and V2(∞) is consistent with V2_0 = 0 trivially.
* `V1_short_coeff` is the leading-order coefficient (1-D skin); subleading tip contributions of
  order (ct/L²)^{3/4} are neglected (estimated small: ratio ~ 2.4·√(ct)/L at ct/L² ≪ 1).
* The total-stress tractions on the faces are transmitted unchanged through the consolidation
  skin (1-D equilibrium), so the outer undrained elastic field is unaltered at O(√t); this and
  the tip subleading estimate are the two structural assumptions behind A.
* Not verified against a full transient numerical solution of the coupled 2-D problem (time
  budget); the uncertain quantity is mainly A (all others are exact elastic/steady limits,
  checked above).
