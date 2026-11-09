# ANSWER.md — First connector failure of a bridged finite slit (plane strain)

## Answers

1. **Remote stress at first connector failure:  Σ = 133.73 MPa**  (133.730 MPa)
2. **Crack advance at each end:  L = 4.086 mm**  (4.0856 mm; total crack half-length c = a + L = 6.4856 mm)

## Method

Quasistatic symmetric growth of the slit from half-length a to c = a + L with
constant closure t = 160 MPa on the zones a < |x| < c (rectilinear bridging law,
"connectors") and intrinsic tip toughness K0 = sqrt(E'Γ0), E' = E/(1−ν²) =
141.059 GPa, K0 = 1.76162 MPa√m. The state on the monotonic branch is the
classical partially-bridged center-crack solution (Tada; Bao & Suo 1992 eqs
(45),(46), generalized to K_tip = K0 ≠ 0; the knowledge unit contains both the
formulas and an independent spectral solver):

- Tip condition (46):  K_tip = √(πc)·[Σ − (2t/π)·arccos(a/c)] = K0
  → **Σ(c) = K0/√(πc) + (2t/π)·arccos(a/c)**  (initiation Σ(a) = K0/√(πa) = 20.29 MPa).
- Opening at the trailing edge x = a (the connector adjacent to the original slit,
  the most-stretched connector since δ decreases toward the tips) (45):
  **δ_t(L) = (4K0/(E'√(πc)))·√(c²−a²) + (8ta/(πE'))·ln(1 + L/a)**.

First connector failure: δ_t(L) = δ_c = 9.00 μm. This 1-D root has a unique
solution (δ_t and Σ both increase monotonically along the branch — no snap-back
before the event):
**L = 4.0856 mm, Σ = 133.73 MPa** (K0-term 12.34 + bridging term 121.39 MPa;
δ_t = 2.109 μm K-part + 6.891 μm bridging part).

## Checks

1. **Two independent COD evaluations agree to 6 digits**: the closed form (45)
   vs. direct adaptive quadrature of the logarithmic-kernel representation
   δ(x) = (2/πE')∫q(s)·ln[(√(c²−x²)+√(c²−s²))/(√(c²−x²)−√(c²−s²))]ds with q = Σ
   on the slit and −t on the zones (solve.py routes A and B: both give
   L = 4.085609 mm, Σ = 133.73029 MPa).
2. **Independent spectral dislocation solver** (knowledge unit `crack_bridging.py`,
   unit-remote + unit-bridging superposition imposing K_tip = K0): N = 200…1600
   gives Σ = 133.66…133.89 MPa, δ_t oscillating around 9.00 μm — consistent
   within its documented ~0.3–0.5% trailing-edge accuracy.
3. **Dugdale/SSB limits** of the formulas: bridging K → −t√(8L/π) (ratio 1.0003 at
   L/a = 4.2e−3); with K0 = 0, δ_t → 8tL/πE' (exact).
4. **K0 = 0 limit reproduces eq. (47)** of the knowledge paper: δ_t = δ_c gives
   Σ = (2t/π)arccos(exp(−π/(8α))) = 131.84 MPa at L = 6.39 mm for α = a/ℓ = 0.302
   (ℓ = E'δ_c/t = 7.93 mm); adding the small intrinsic toughness (g = Γ0/tδ_c =
   0.0153) shifts the event to smaller L, higher Σ, as expected.
5. **Physics sanity**: Σ = 133.7 MPa lies between initiation (20.3 MPa) and the
   large-L bridging ceiling Σ → t = 160 MPa; K_applied = Σ√(πc) = 19.07 MPa√m is
   carried mostly by the bridged zone (K_bridge = 17.31 MPa√m), tip keeps
   K_tip = K0 = 1.76 MPa√m exactly.

## Code

`solve.py` (this folder) — constants, closed-form route, quadrature route,
spectral-solver route, validations 1–5, root find via `brentq`.
