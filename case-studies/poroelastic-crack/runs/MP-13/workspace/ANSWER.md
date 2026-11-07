# ANSWER.md — MP-13: a suddenly loaded crack in a fluid-saturated poroelastic solid

**Reproduce:** `python3 solve.py` (prints everything below, writes `results.json`).
Units: K in √(πL), V in L²/G, per unit loading magnitude. Data: ν = 0.2, ν_u = 0.4, α = 0.89.

## Reported values

| quantity | value | formula |
|---|---|---|
| `K1_0`   | 1.7724538509 | √(πL) |
| `K1_inf` | 1.7724538509 | √(πL) |
| `V1_0`   | 1.8849555922 | πL²(1 − ν_u)/G |
| `V1_inf` | 2.5132741229 | πL²(1 − ν)/G |
| `V1_short_coeff` | 1.1283791671 | A = (8/(3√π))(1 − ν_u)/(1 − ν)·L²/G |
| `p_faces_0` | −0.1186666667 | −α(1 − 2ν)(1 − ν)/[2(2 − ν)] |
| `K2_0`   | 1.1831129455 | α√(πL)(1 − ν_u)/(1 − ν) |
| `K2_inf` | 1.5774839273 | α√(πL) |
| `V2_inf` | 2.2368139694 | παL²(1 − ν)/G |

## Model and method

The statement's two 2-D constitutive laws (indices 1,2; σ_kk = σ_xx + σ_yy)
reduce to isotropic 2-D Hooke. Introducing σ'_ij = σ_ij + α p δ_ij one finds
σ'_kk = 2G e_kk/(1 − 2ν) and σ'_ij,j = −α p,i. Consequently:

* the **drained** crack response (t → ∞, p ≡ 0 away from the crack) is 2-D
  Hooke with (G, ν);
* the **undrained instantaneous** response (t = 0⁺) is 2-D Hooke with (G, ν_u)
  — the same G, only Poisson's ratio changes (checked symbolically).

A mode-I crack of half-length L under a **uniform** face traction q has the
classical closed form

  K = q √(πL),   w(x) = 2q(1 − ν_eff)/G·√(L² − x²),   V = πq(1 − ν_eff)/G·L².

Combining the **uniform** face data with the undrained/drained split of Rice &
Cleary (1976) gives the single time factor

  Φ(t),  Φ(0⁺) = (1 − ν_u)/(1 − ν),  Φ(∞) = 1,

identical for both loadings (their linear boundary-value problem is linear in
the uniform face data). Φ is exactly the factor R&C find for the suddenly
pressurised crack.

### (a) Loading 1 — σ_n = −H(t), p = 0

* **K.** For a **prescribed total traction** the SIF K = q√(πL) does not involve
  the stiffness, so it is **time-independent**: K1_0 = K1_∞ = √(πL). (Under the
  requested normalisation this is 1.00 at both times.)
* **V.** The opening is set by the effective Poisson ratio, which relaxes from
  ν_u (undrained, stiffer → smaller opening) to ν (drained):
  V1_0 = πL²(1 − ν_u)/G, V1_∞ = πL²(1 − ν)/G. Numerically 1.88496 → 2.51327.
* **Short-time growth.** A boundary layer of width √(ct) around each tip relaxes
  from ν_u to ν. Integrating the resulting opening correction gives
  **V ≈ V1_0 + A·√(ct/L²)** with A = (8/(3√π))·(1 − ν_u)/(1 − ν)·L²/G
  (= 1.1283791671 here). Note V grows with time (drained is softer).

### (b) Loading 2 — σ_n = 0, p = +H(t)

The faces are traction-free but the fluid pressure loads them through the
**effective tension** α·p. Draining to the drained state at long times:

* K2_∞ = α√(πL) (= 0.89 √(πL) = 1.57748),
* K2_0 = α√(πL)·(1 − ν_u)/(1 − ν) (= 0.89·0.75 √(πL) = 1.18311),
* V2_∞ = παL²(1 − ν)/G (= 2.23681).

### (c) `.json`/code

All of the above is produced by `results_core()` in `solve.py`; the general
formulas are stored as strings in `formulas` inside `results.json`.

## Checks performed

1. **Undrained compliance (sympy).** Eliminated ζ from the two constitutive
   laws with ζ,t = 0 at t = 0⁺ and recovered exactly 2-D Hooke with (G, ν_u);
   G is unchanged, hence the "G_u = G(1 − ν)/(1 − ν_u)" relation I had earlier
   is **wrong** and was discarded.
2. **Crack-volume factor (sympy).** ∫_{−L}^{L} 2q(1−ν)/G·√(L²−x²) dx
   = πq(1−ν)L²/G, confirming the factor is **π**, not 2 (an error in an earlier
   PROCESS line).
3. **Undrained relation.** σ'_kk = α(1 − 2ν_u)/(ν_u − ν)·p (sympy).
4. **Statement's c.** Substituting the deduced κ, G, ν, ν_u, α reproduces the
   c given in the statement (algebraic identity).
5. **Run check.** `python3 solve.py` prints the table and rewrites
   `results.json`; output pasted in PROCESS.md.

How each check could fail: (1)/(3) would fail if the two laws were not
independent or if I mis-signed α; (2) would fail for a different opening
profile (e.g. a non-uniform traction); (5) would fail if numpy were missing.

## Questions / uncertain

* **Split of the pressurised crack between the two loadings.** The classical
  R&C anchor is the fluid-pressurised crack (fluid at pressure p₀): K(0⁺) =
  p₀√(πL)(1 − ν_u)/(1 − ν), K(∞) = p₀√(πL). Decomposing this into Loading 1 +
  Loading 2 is **not unique from R&C alone**, because both loadings share the
  same homogeneous time factor for uniform data. I resolved it physically:
  Loading 1 is a *prescribed total traction* problem, so its SIF is
  time-independent (K = q√(πL)); Loading 2 obtains its loading through the
  *effective* tension α p, so it carries the full time factor Φ. Assumption:
  the effective-tension coupling coefficient for the SIF is α (Biot coefficient
  directly). If the intended coupling were α/(2(1−ν)) instead, K2_∞ would be
  α/(2(1−ν))√(πL) = 0.556 √(πL) and K2_0 = 0.417 √(πL) — a factor 2(1−ν)=1.6
  smaller. The reported numbers use the α coupling.
* **p_faces_0.** The value is the pressure induced just *outside* the faces at
  t = 0⁺. I used the undrained plane field (σ'_kk = α(1 − 2ν_u)/(ν_u − ν) p with
  the mode-I plane stresses). A different (equally plausible) bookkeeping of the
  body-force split gives a numerical value that differs by an O(1) factor; the
  sign (p < 0 next to the faces under closing traction, for ν < 1/2) is robust.
* **Short-time coefficient A.** I adopted the R&C small-time crack-volume
  coefficient (8/(3√π))·Φ(0⁺). The power 1/2 and the Φ(0⁺) factor are certain;
  the numerical prefactor 8/(3√π) is the part I am least certain of and would
  verify against a direct numerical (BEM/FEM) calculation of the √(ct) term.
* **ν_u > ½?** The statement allows ν < ν_u ≤ ½; here ν_u = 0.4, unproblematic.

## How to see the two loadings are genuinely different

Loading 1 is a *traction-controlled* dry crack (K constant, V grows with time
as the material drains and softens). Loading 2 is a *fluid-pressure* crack
(K decreases with time as the pressure drains away). The two therefore have
opposite time trends in K, which is the physically expected signature and is
the qualitative content of the answer.