# ANSWER — a suddenly pressurized crack in a fluid-saturated porous solid

Reproduce with: **`python3 solve.py`** (prints every value in `results.json`
and rewrites it).

## 0. Summary of the answers (part (c): ν = 0.2, ν_u = 0.4, α = 0.89)

K in units of √(πL), V in units of L²/G, per unit loading magnitude.

| quantity | symbol/formula | value |
|---|---|---|
| `K1_0` | `1·√(πL)` | 1.0000000000 |
| `K1_inf` | `1·√(πL)` | 1.0000000000 |
| `V1_0` | `2π(1−ν_u)L²/G` | 3.7699111843 |
| `V1_inf` | `2π(1−ν)L²/G` | 5.0265482457 |
| `V1_short_coeff` | `(V1_inf−V1_0)·2/√π` | 1.4179630807 |
| `p_faces_0` | `2(ν_u−ν)/[α(1−2ν)(1−ν_u)]` | 1.2484394507 |
| `K2_0` | `0` | 0.0 |
| `K2_inf` | `0` | 0.0 |
| `V2_inf` | `2π(1−ν)αL²/G = α·V1_inf` | 4.4736279387 |

## 1. Model and method

**Model.** Biot linear poroelasticity exactly as given in the problem:
the constitutive law (1) `2G e_ij = σ_ij − ν σ_kk δ_ij + α(1−2ν) δ_ij p`,
the ζ relation (2), equilibrium, Darcy's law and mass balance. Plane strain,
infinite medium, crack −L ≤ x ≤ L on y = 0, step loadings for t > 0.

**Method.**
1. Identify the *effective stress* hidden in law (1) so that the crack
   boundary-value problem reduces to an ordinary (drained) elastic crack
   problem, and identify the undrained limit `ζ = 0`.
2. Use the classical mode-I elastic crack solution (uniform traction on the
   faces) to get **K** and the **opening profile / volume**.
3. For the time dependence of loading 1, use the **1-D pore-pressure
   diffusion** of the undrained jump next to the drained crack face.
4. Verify every coefficient with independent numerics (quadrature of the
   opening, a finite-difference solve of the diffusion problem, and symbolic
   algebra for the algebraic identities and the surface pore pressure).

## 2. Derivation in brief

### 2.1 Effective stress and the drained / undrained limits

Seek `σ^eff_ij = σ_ij + β p δ_ij` (2-D indices). Then
`σ^eff_ij − ν σ^eff_kk δ_ij = σ_ij − ν σ_kk δ_ij + β(1−2ν) p δ_ij`.
Comparing with law (1) gives `β = α`, hence

> **`σ^eff_ij = σ_ij + α p δ_ij`** (standard Biot effective stress), and law (1)
> *is* Hooke's law in `σ^eff` with the **drained** ν.

(Checked symbolically: `verify_effective_stress.py`, residual = 0.)

Setting `ζ = 0` (undrained) in (2) gives
`p = −(ν_u−ν)/[α(1−2ν)] · σ_kk`; substituting into (1) collapses it to
`2G e_ij = σ_ij − ν_u σ_kk δ_ij`, i.e. **Hooke with ν_u**
(`derive.py`, `verify_V1_0.py`, residual = 0). So:

* `t = 0⁺` (no fluid has moved, ζ = 0): elastic response with **(G, ν_u)**.
* `t → ∞` (faces drained, p ≡ 0): elastic response with **(G, ν)**.

### 2.2 Stress intensity factors

For a crack in an infinite body loaded by a **prescribed traction** on its
faces, `K_I` is **independent of the elastic constants** (statically
determinate): `K = P √(πL)` for a uniform traction `P`.

* Loading 1: total face traction `σ_n = −1` at all times ⇒
  **`K1_0 = K1_inf = 1·√(πL)`** (the undrained/drained distinction does not
  change K).
* Loading 2: total face traction `σ_n = 0` ⇒ no mode-I driving force ⇒
  **`K2_0 = K2_inf = 0`**.

### 2.3 Crack volume (loading 1)

The uniform-traction mode-I opening for face traction `P` is
`w(x) = 4(1−ν)P/G · √(L²−x²)` (total displacement jump; verified by exact
mode-I stress field and by quadrature, `check_elastic.py`: the volume
quadrature matched `2π(1−ν)PL²/G` to 4e−10). Integrating,

> `V = ∫_{−L}^{L} w dx = 2π(1−ν) P L² / G`.

Therefore
* **`V1_0 = 2π(1−ν_u) L²/G`** (undrained, t = 0⁺),
* **`V1_inf = 2π(1−ν) L²/G`** (drained, t → ∞).

### 2.4 Small-time growth of V (loading 1)

Just after loading, the crack faces are held at `p = 0` while the solid just
next to them has the undrained jump `p0` (§2.5). The pressure deficit
`u = p0 − p` diffuses normal to the face (`∂u/∂t = c ∂²u/∂y²`,
`u(0,t)=p0`, `u(y,0)=0`) with the exact solution `u(y,t)=p0·erfc(y/√(4ct))`,
i.e. **`p(y,t) = p0·erf(y/√(4ct))`** (`diffusion.py`, verified against an FD
solve). The crack volume relaxes linearly in the induced pressure, from the
undrained value at t = 0⁺ to the drained value as t → ∞:

`V(t) = V1_inf − (V1_inf − V1_0)·S(t)`, `S(t) = (1/L)∫_0^L p(y,t)/p0 dy`.

For `ct/L² ≪ 1`, `S(t) = 1 − (2/√π)√(ct)/L + …`, hence

> **`V(t) = V1_0 + A √(ct/L²)`**, **`A = (V1_inf − V1_0)·2/√π = 4√π(ν_u−ν)L²/G`**.

Numerically (`check_A.py`): the slope of `V − V1_0` vs `√(ct/L²)` converged to
`1.4179645…` against the closed form `1.4179630807` as `ct/L² → 10⁻⁸`.

### 2.5 Induced pore pressure next to the faces at t = 0⁺ (loading 1)

The exact mode-I stress field on the crack faces gives
`σ_yy = σ_xx = −P` (`σ_xy = 0`), i.e. the in-plane stress at a face point is
hydrostatic of magnitude P; verified numerically for `|x| < L`
(`mode1_stress.py`: `σ_xx = σ_yy = −0.99999…` up to the tip). At `t = 0⁺`
the medium is undrained (`ζ = 0`) and plane strain (`e_zz = 0`). With
`σ_xx = σ_yy = −1`, `ζ = 0`, `e_zz = 0`, solving (1)+(2) for `(σ_zz, p)`
(`derive.py`) gives

> **`p_faces_0 = 2(ν_u − ν) / [α (1−2ν)(1−ν_u)]`** = 1.2484394507.

(The denominator factors as `(2ν−1)(ν_u−1)`; the sign is positive, i.e.
compressive pore pressure, consistent with the compressive face load.)

### 2.6 Loading 2

Faces are traction-free with pore pressure `p = 1`. By the effective-stress
result (§2.1) the crack opens under the **effective** face tension
`σ^eff_n = σ_n + α p = α`. Hence, in the long-time (drained) limit,

> **`V2_inf = 2π(1−ν) α L² / G = α V1_inf`** = 4.4736279387.

`K2_0 = K2_inf = 0` since the total face traction is zero throughout.

## 3. Checks performed (what was compared, output pasted, how it could fail)

1. **Elastic opening profile & volume** (`check_elastic.py`): compared the
   closed-form `w(x)=4(1−ν)P/G √(L²−x²)` against a fine-grid quadrature.
   Output:
   ```
   V numeric = 5.0265482438620985
   V analytic= 5.026548245743669
   rel err   = 3.7432654836073345e-10
   [tip]  f=1e-08  w/sqrt(2L(L-x)) = 3.1999999992   4(1-nu)P/G=3.2000000000
   ```
   *Could fail if* the opening prefactor were `2` instead of `4`; the
   quadrature would then be off by a factor 2.
2. **Effective stress / undrained reduction** (`verify_effective_stress.py`,
   `derive.py`, `verify_V1_0.py`): symbolic residuals exactly `0`:
   ```
   2-D effective-stress residual (should be 0): 0
   undrained reduction residual (should be 0): 0
   ```
   *Could fail if* the `δ_ij` trace were 3-D instead of 2-D; then the residual
   would be `−ανp ≠ 0` (exactly what happened on my first, wrong attempt —
   see PROCESS.md — and fixed here).
3. **Face stress state** (`mode1_stress.py`): exact Westergaard mode-I field:
   ```
   x=0.50 y=1e-09: sxx=-1.000000  syy=-1.000000  sxy= 0.00e+00
   ```
   *Could fail if* `σ_xx` at the face were load/geometry dependent; the
   `p_faces_0` formula would change.
4. **Diffusion profile** (`diffusion.py`): FD solve vs exact `erf`:
   ```
   t=0.500
      y=0.0000  FD=1.000000  exact=1.000000 diff= 0.00e+00
      y=0.5942  FD=0.552544  exact=0.548506 diff= 4.04e-03
   ```
   Boundary value exact; interior within the (coarse) FD resolution; the
   discrepancy shrinks with finer `dy`. *Could fail if* the face BC were
   undrained rather than drained (then no decay; `V` would not relax to
   `V1_inf`).
5. **Small-time coefficient** (`check_A.py`):
   ```
   ct/L^2=1e-07  S=0.9996431751   slope=(V-V0)/sqrt(ct/L^2)=1.4179632284
   (R1) closed-form  A = (Vinf - V0)*2/sqrt(pi) = 1.4179630807
   ```
   *Could fail if* the diffusion length scale were not `L` (an O(1) prefactor
   change); see §4 item 2.
6. **V1_0** (`verify_V1_0.py`):
   ```
   V1_0 numeric = 3.769911182896573  analytic = 3.7699111843077517  rel err = 3.743267839573765e-10
   ```
7. **Whole-solution reproducibility** (`python3 solve.py`): prints
   `K1_0=1.0, K1_inf=1.0, V1_0=3.7699111843077517, ...` matching the table in
   §0 and `results.json`.

## 4. Questions / assumptions / uncertainties

1. **Convention for "crack volume per unit thickness".** The problem defines
   it as "the integral of the opening over the crack". I used the total
   displacement jump `w = 4(1−ν)P/G√(L²−x²)`, giving `V = 2π(1−ν)PL²/G`.
   *Assumption.* If "opening" meant the one-sided half-opening, all V's would
   be halved; ratios such as `V2_inf/V1_inf = α` are unaffected.
2. **Influence depth in the small-time coefficient A.** I took the
   opening-influence depth to be `L` (the only length in the boundary-layer
   problem), yielding the classical `2/√π` coefficient. A more elaborate
   self-consistent poroelastic weight function could change the O(1)
   prefactor; the `√(ct/L²)` *exponent* and the proportionality
   `A ∝ (ν_u−ν)L²/G` are robust. *Assumption: coefficient 2/√π.*
3. **`p_faces_0` "immediately next to the crack faces".** Interpreted as the
   undrained value at the face approached from the solid side (no fluid
   motion yet). Value `1.248` exceeds the applied traction `1`, a consequence
   of the low drained ν combined with poroelastic stiffening. *Assumption.*
4. **K for loading 2.** Interpreted literally: total face traction is zero, so
   `K = 0`. If a "hydraulic" SIF driven by the effective pressure `αp` were
   intended, it would be `K = α√(πL)` in this convention; the statement says
   "total stress", so `K2 = 0`.
5. **`V1_short_coeff` sign.** `ν_u > ν` ⇒ `A > 0` ⇒ the crack opens further
   with time, as expected for a drained face under fixed traction.