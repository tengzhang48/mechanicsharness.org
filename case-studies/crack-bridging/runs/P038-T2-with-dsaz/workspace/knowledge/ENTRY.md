# ENTRY.md — Crack-bridging concepts (Bao & Suo 1992, Appl. Mech. Rev. 45:355)

## Problem

Predict the mechanical properties of a brittle-matrix composite (or monolithic
solid) from a **bridging law** σ(δ) that relates the bridging (closure) stress
across a matrix crack to the crack-face separation. The law has a strength
scale σ₀ and a limiting separation δ₀, which combine with the plane-strain
modulus E′ to define the single material length

    ℓ = δ₀ E′ / σ₀        (ranges from 10⁻⁹ m for atomic bonds to ~1 m for
                           cross-over fibers; Table 1)

Small-scale bridging (SSB): damage zone L₀ ≪ all crack/component lengths; the
exterior problem reduces to a semi-infinite crack in a K-field; toughness Γ and
R-curves are well defined. Large-scale bridging (LSB): damage zone comparable
to a component length a; toughness concepts fail; strength must be computed
from the coupled component + bridging-law problem; the controlling group is

    α = a / ℓ.

## Assumptions and validity range

- Bridging law σ = σ₀χ(δ/δ₀), eq (1): separates the *shape* χ from the *scale*
  (σ₀, δ₀); E′ = E (plane stress) or E/(1−ν²) (plane strain). [p. 355, 357]
- Linear spring law σ = σ₀ + sδ for δ<δ₀, 0 beyond: s>0 hardening, s<0
  softening; β = sL/E′ is the only dimensionless parameter of the SSB problem.
  [p. 359]
- Steady (plateau) state of an R-curve: separation at the trailing edge of the
  bridging zone reaches δ₀; tip condition K_tip = K₀ (matrix toughness). [p. 359]
- LSB strength: crack of half-length c = a+L in an infinite panel, uniform
  closure on a<|x|<c, failure at δ_t = δ₀ with K_tip = Γ₀ = 0; valid when
  σ(L) increases monotonically (fails for softening laws with interior
  maximum — then maximize σ over L with δ_t ≤ δ₀). [p. 361–362]
- All SSB solutions below are plane-strain/plane-stress agnostic (only E′
  enters). Inherent assumption: no fiber failure before δ₀ (pull-out law (2)
  is a conservative design limit).

## Main results

### Material length and energies (Section 1)
- Γ = σ₀δ₀ and L₀ = (π/8)ℓ for the rectilinear law, eqs (10),(11).
  [reproduced: exact identity, eq (19) coefficient c1 = √(8/π) = 1.59577 from
  the zone solver to 2e-5; L₀ = π/8 = 0.39270]
- Table 1 values σ₀δ₀ = (1, 10⁴, 10³) J/m² and ℓ = (10⁻⁹, 10⁻³, 1) m for
  (atomic bond, fiber pull-out, fiber cross-over) at E = 10¹¹ N/m².
  [reproduced: all printed products follow exactly from the printed σ₀, δ₀;
  verified against the page image — the table is internally consistent, no
  erratum. An earlier draft misread δ₀ = 10⁻⁹ and invented a 10× energy
  erratum; withdrawn.]
- Fiber pull-out, eqs (2)–(4): σ₀ = fS_b, δ₀ = (1−f)²S_b²R/(2Eτ), hence
  Γ ~ f(1−f)²S_b³R/(2Eτ), L₀ ~ π(1−f)²S_bR/(16 fτ). [reproduced: follows
  algebraically from (3),(4) + (10),(11) — used as a consistency check]

### Small-scale bridging solutions (Section 2)
- Rectilinear law: Γ = Γ₀ + σ₀δ₀ (18); R-curve K = √(8/π)σ₀√L + K_tip (19);
  L₀(Γ₀) = (π/8)ℓ[√(1+g) + √g]⁻², g = Γ₀/(σ₀δ₀) (20).
  [reproduced: (18),(19) exact; (20) checked to ~5e-3 by a numerical
  steady-state construction that solves the zone integral equation at each
  trial L (no closed-form coefficients used); the printed −2 exponent is
  forced — a +2 bracket would be off by a factor
  [(√(1+g)+√g)/(√(1+g)−√g)]², and it is the form that makes L₀ decrease
  as Γ₀ grows]
- Linear spring law: K = c₁σ₀√L + c₂K₀ (23), δ_t = c₃σ₀L/E′ + c₄K₀√L/E′ (24),
  with c₂ = √(1+2βc₁²), c₃ = (c₂−1)/β, c₄ = 2c₁ (26).
  [reproduced: all three identities hold to 2e-3 for β ∈ (−0.37, +1.6);
  Table 2 (hardening): c₁ matched to ≤ 0.005, c₂ (from (26)) to ≤ 0.009 at all
  ten tabulated 8β/π (measured maxima 0.0040 / 0.0083)]
- Softening triangle σ = σ₀(1−δ/δ₀): plateau Γ = ½σ₀δ₀ (27), L₀ = 0.366 ℓ
  (28); R-curve 2𝒢/(σ₀δ₀) = 2c₁²λ, λ = L/ℓ; Table 3.
  [reproduced: all ten rows of Table 3 to ≤ 0.003 (c₁) / 0.002 (2𝒢/σ₀δ₀),
  measured maxima 0.0020 / 0.0008; plateau λ* = 0.3660 ± 2e-3,
  𝒢_plateau/(σ₀δ₀) = 0.5000 ± 5e-3, c₁(λ*) = 1.1688 vs 1.169 printed]
- Power-law bridging σ/σ₀ = (δ/δ₀)^m: scaling forms (30)–(32) with
  λ = L/(δ₀E′/σ₀). [taken: dimensionless functions Δ, k tabulated elsewhere
  (Gu et al. 1992), not in this paper — not reproduced]

### Large-scale bridging consequences (Section 3)
- Griffith/SSB strength of a monolithic ceramic, eq (44):
  S/σ₀ = (πα)^{-1/2}, α = a/ℓ. [reproduced: closed form taken; consistent
  with LSB curve at large α]
- LSB strength with rectilinear bridging, eqs (45),(46) (Tada) and (47):
  (8/π)α ln[sec(πS/2σ₀)] = 1, i.e. S/σ₀ = (2/π)arccos(e^(−π/(8α))).
  [reproduced two ways: (47) derived exactly from (45)+(46) since the
  bracket of (45) vanishes at K_tip = 0; and the spectral dislocation solver
  agrees with (47) to 0.06–0.65 % at α = 0.2…1.4]
- Strength curves from bridging-law shape, Fig. 11: with the same σ₀, δ₀ the
  softening triangle gives a lower strength curve than the rectilinear law,
  e.g. at α = 0.2: 0.65 vs 0.91; α = 1: 0.36 vs 0.53 (read from figure;
  softening values reproduced by the spectral solver within ~0.02).
- Notch sensitivity, eqs (48),(49) + Fig. 13: notch-brittle σ̄_max = σ_b/C
  (C = 3 for a small hole, a/w→0, hence 1/3 asymptote), notch-ductile
  net-section σ̄_max = σ_y(1−a/w). [reproduced: limiting curves exact; the
  intermediate bridging curves of Fig. 14 require FE computations (Ho & Suo
  1992) — not reproduced]
- Splitting beam, eqs (51)–(54): J = 12M²/E′h³ (independent of the bridging
  law — the special LSB specimen), M_max from δ_t = δ₀, L₀ = a₁²(E′δ₀/3σ₀)h
  with a₁(L/h) of order unity. [taken, not reproduced numerically]
- Measuring bridging laws via J = J_tip + ∫σ dδ (56),(57): valid under SSB;
  under LSB J is not related to applied load through handbook K — except for
  special geometries like the splitting beam (58). [conceptual]

## Conditions and limits

- SSB ⇔ α ≫ 1 (eq (14)); LSB when α ≲ 1. The same group α takes different
  meanings: unbridged crack size (toughness), flaw size (strength), hole radius
  (notch). [p. 358, 361]
- Toughness/K/R-curves are meaningful only under SSB; "effective" R-curves
  converted under LSB are geometry-dependent and misleading. [p. 358]
- The strength procedure (δ_t = δ₀ at K_tip = 0) presumes σ(L) monotonic;
  for softening laws or Γ₀ ≠ 0 the maximum may occur earlier (Carpinteri 1990,
  Bao & Zok 1992) — use the constrained maximum. [p. 361]
- Eq (26) holds exactly for the linear spring law only; different shapes
  modify prefactors within order unity (Figs 8, 11 bounds). [p. 360–362]
- Table 1 is illustrative (order-of-magnitude) values; cross-over fiber data
  inferred from CAS/SiC experiments.

## Relations to other results

- Eqs (10),(11),(19),(20) are the classical Dugdale–Bilby–Cottrell–Swinden
  strip results; the zone solver here reproduces the Dugdale identities
  K_tip = −σ₀√(8L/π) and δ_t = 8σ₀L/(πE′) exactly.
- Eq (44) is the Griffith result recast: with Γ = σ₀δ₀, K²/E′ = πaσ²/E′.
- Eq (47) is the fully-bridged limit of the partially-bridged crack solution
  (45),(46); it matches the small-α asymptote S/σ₀ → 1 (strength bounded by
  the closure strength) and approaches eq (44) from below as α → ∞ without
  ever crossing (relative gap ≈ 9 % at α = 0.7, ≈ 3 % at α = 2, < 1 % for
  α ≳ 5; Fig. 10 shows an approach, not an intersection).
- The normalization idea (shape vs scale, single length ℓ) parallels
  Needleman (1990) cohesive laws; the LSB dimensionless group α is the
  analogue of the ductile-fracture group a/(J/σ₀²).
- The softening-triangle coefficients (Table 3) come from Bao & Hui (1990);
  the hardening ones (Table 2) from Budiansky, Amazigo & Evans (1988); eqs
  (45),(46) from Tada, Paris & Irwin (1985).

## Errata / inconsistencies found

**None in the paper.** After checking the page images at 300 dpi, every
quantitative item transcribed in this unit is internally consistent.

1. **Table 1 (withdrawn transcription error, corrected):** an earlier draft
   of this unit read the atomic-bond δ₀ as 10⁻⁹ m (and ℓ as 10⁻⁸ m) and
   flagged the printed σ₀δ₀ = 1 J/m² as an erratum. The page image
   (`pages/page-02.png`) clearly prints δ₀ = 10⁻¹⁰ m, which makes
   σ₀δ₀ = 1 J/m² and δ₀E/σ₀ = 10⁻⁹ m exactly. The paper is consistent;
   the error was in our transcription. [resolved]
2. **Eq (20) misprint risk** (extraction note, resolved): the bracket
   carries an exponent −2,
   L₀ = (π/8)ℓ[√(1+g) + √g]⁻². Reading it as +² contradicts eqs (18)+(19)
   and the paper's own statement "the tougher the matrix, the smaller the
   bridging zone". [verified numerically against the zone-solver route]
3. Table 2, row 8β/π = 1: c₂ printed as "1.91" (two decimals); value
   consistent with eq (26). [cosmetic]
4. Text eq for (13) in the OCR'd text is garbled ("δ_b"); the printed
   equation reads S_b (fiber bundle strength). [extraction note]

## Using this unit

**Questions this unit answers**
- What are the SSB (small-scale bridging) toughness and R-curve for a given
  bridging law? (rectilinear: eqs 18–20; linear spring: eqs 23–26 + Table 2;
  softening triangle: eqs 27–28 + Table 3)
- What is the steady-state bridging zone size L₀ with/without matrix
  toughness Γ₀? (eqs 11, 20, 28)
- What strength does a bridged crack of size a have when the bridging zone
  is comparable to a (LSB)? (eqs 44, 47; Fig. 11 shapes)
- How notch-sensitive is the solid? (eqs 48, 49 limiting curves)
- What material length ℓ and energy σ₀δ₀ correspond to a bridging mechanism?
  (Table 1, eqs 10–13)

**Input conventions and units**
- Everything is normalized by the bridging scales: set σ₀ = δ₀ = E′ = 1 in
  the solvers; dimensionless groups are β = sL/E′ (spring stiffness),
  λ = Lσ₀/(E′δ₀), g = Γ₀/(σ₀δ₀), α = a/ℓ with ℓ = E′δ₀/σ₀.
- Physical units only enter `table1()` (SI: N/m², m, J/m², E = 10¹¹ N/m²).
- Plane stress vs plane strain enters only through E′ (E vs E/(1−ν²)).

**Functions to call (all in `crack_bridging.py`; zone solver in `zone.py`)**
- `table1()` → dict of Table 1 inputs and printed products/lengths.
- `c1_rectilinear(L=0.01)` → eq (19) coefficient √(8/π) from the zone solver.
- `L0_rectilinear(g)` → eq (20) closed form, units of ℓ.
- `steady_state_rectilinear(g, N=48)` → same quantity, numerical route
  (solves the zone integral equation; independent of the closed form).
- `hardening_coefficients(beta)` / `reproduce_table2()` → eqs (23)–(26)
  coefficients c₁…c₄ / full Table 2 reproduction.
- `softening_rcurve_point(lam)` / `reproduce_table3()` → Table 3 R-curve
  points (c₁, 2𝒢/σ₀δ₀, δ_t/δ₀).
- `softening_plateau()` → eqs (27),(28): λ*, c₁(λ*), 𝒢_plateau/(σ₀δ₀).
- `rcurve_rectilinear_x(x)`, `rcurve_softening_points()` → Fig. 8 curves.
- `strength_SSB(alpha)`, `strength_LSB(alpha)` → eqs (44), (47) closed forms.
- `strength_numerical(alpha, s=0.0 or -1.0, N=700, nL=30)` → spectral-solver
  strength (s=0 rectilinear, s=−1 softening triangle); maximizes σ(L) s.t.
  δ_t ≤ δ₀, so it also handles softening laws with an interior maximum.
- `notch_brittle_limit(0.0)` → 1/3 (eq 48, a/w→0 only);
  `notch_ductile_limit(a/w)` → 1−a/w (eq 49).

**Where the model holds — and where it does not (physical reasons)**
- SSB results (toughness, R-curves, Tables 2–3) require the bridging zone
  ≪ crack/component size (α ≫ 1). The semi-infinite solvers assume this;
  for α ≲ 1 the exterior boundary couples to the zone and K-concepts fail.
- LSB strength (eqs 44, 47) is for a center crack in an infinite panel,
  rectilinear bridging, Γ₀ = 0, and σ(L) increasing monotonically; for
  softening laws use the constrained maximum implemented in
  `strength_numerical` (the paper's own caveat, p. 361).
- Eq (26) identities hold exactly only for the linear spring law; other
  shapes modify prefactors by order-unity factors.
- The solvers assume no fiber failure before δ₀ (pull-out law (2) is a
  conservative design limit) and a single bridging scale (no law-shape
  transition inside the zone).
- Notch: only the two limiting curves are analytic; intermediate curves need
  FE (Ho & Suo 1992). Splitting-beam numbers a₁(L/h) are not in the paper.

**Provenance of the main results**
- Taken from the paper (closed forms): eqs (18)–(20), (23)–(28), (44),
  (47)–(49), Tables 1–3 transcriptions.
- Re-derived here: eq (47) from (45)+(46) (bracket vanishes at K_tip = 0);
  eqs (12),(13) from (3),(4)+(10),(11).
- Implemented (numerical): spectral dislocation solver (finite center crack)
  and semi-infinite zone solver; `strength_numerical`; the numerical route
  to eq (20).
- Independently checked: zone solver vs Dugdale identities (K_tip = −σ₀√(8L/π),
  δ_t = 8σ₀L/(πE′), exact); spectral solver vs uniform-pressure and
  remote-loading closed forms (1e-12); eq (20) numerical route vs closed
  form (5e-3); eq (47) vs spectral solver (0.06–0.65 %); Tables 2–3 vs the
  printed digits (≤0.004/0.009 and ≤0.002/0.001); Table 1 arithmetic.
