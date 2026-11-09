# REPORT.md — Reproduction of Bao & Suo (1992) "Remarks on Crack-Bridging Concepts"

## What was reproduced (28/28 tests pass, `python3 -m pytest`)

| Item | Paper | Ours | Status |
|---|---|---|---|
| Eq (19) coefficient c₁ = √(8/π) | 1.59577 | 1.59577 (±2e-5) | reproduced |
| Eqs (10),(11) Γ, L₀ (rectilinear) | Γ=σ₀δ₀, L₀=π/8·ℓ = 0.3927ℓ | exact identities | reproduced |
| Eq (20) L₀(Γ₀), g=0…3 | (π/8)ℓ[√(1+g)+√g]⁻² | ≤ 5e-3 (numerical route, not using the closed-form coefficients) | reproduced |
| Eq (26) identities c₂²=1+2βc₁², c₃=(c₂−1)/β, c₄=2c₁ | — | ≤ 2e-3 for β ∈ (−0.37, 1.6) | reproduced |
| Table 2 (hardening c₁, c₂), 10 rows | 1.670…2.414 | max abs error 0.0040 / 0.0083 (tests enforce 0.005 / 0.009) | reproduced |
| Eqs (27),(28) softening plateau | Γ=½σ₀δ₀, L₀=0.366ℓ | 0.5000±5e-3, 0.3660±2e-3 | reproduced |
| Table 3 (softening c₁, 2𝒢/σ₀δ₀), 10 rows | 1.564…1.169, 0.192…1.000 | max abs error 0.0020 / 0.0008 (tests enforce 0.003 / 0.002) | reproduced |
| Fig. 8 kinks | 0.393, 0.732 | π/8 = 0.3927; 2×0.366 = 0.732 | reproduced |
| Eq (44) SSB strength | (πα)^{-1/2} | closed form + consistency | reproduced (taken) |
| Eq (47) LSB strength | implicit curve | closed form derived from (45)+(46); spectral solver agrees 0.06–0.65 % | reproduced |
| Fig. 10 relation | LSB approaches SSB from below | gap 9.0 % (α=0.7) → 0.65 % (α=10); never crosses | reproduced |
| Fig. 10/11 strength curves | points read from figures | rectilinear: matches; softening solver 0.642/0.446/0.361/0.313 vs read 0.65/0.44/0.36/0.33 | reproduced |
| Eqs (48),(49) + Fig. 13 limits | 1/3 at a/w→0; 1−a/w | exact | reproduced |
| Table 1 | σ₀δ₀, ℓ for 3 mechanisms | all printed products exact from printed inputs; verified on page image | reproduced (no erratum) |

Two independent solution methods were built and cross-validated:
1. **Spectral dislocation solver** (`crack_bridging.py`): center crack |x|<c,
   Chebyshev expansion of the dislocation density (kernel exact to machine
   precision for smooth tractions; validated against uniform-pressure and
   remote-loading closed forms to 1e-12).
2. **Semi-infinite zone solver** (`zone.py`): integral equation on the zone
   with the logarithmic kernel (4/πE′)ln|(√x+√ξ)/(√x−√ξ)| and the K-kernel
   √(2/πξ), calibrated against solver 1 (g₀ = 4.015 ± 0.015, residuals 4e-5;
   Dugdale identities exact). Used for all SSB tables (avoids the traction-jump
   singularity at the zone trailing edge that limits solver 1 to ~0.5 % there).

## What was NOT reproduced and why

- **Power-law bridging functions Δ, k (eqs 30–32)**: tabulated elsewhere
  (Gu et al. 1992), not given in this paper. Only the scaling forms were taken.
- **Fig. 14 notch ductile-to-brittle transition curves**: computed with finite
  elements by Ho & Suo (1992); only the two limiting curves (eqs 48, 49) are
  reproducible analytically.
- **Splitting-beam numbers (eqs 51–54)**: the J-integral identities are taken;
  the coefficient a₁(L/h) comes from Suo et al. (1992a) and is not given here.
- **Eqs (12),(13) microstructural scalings**: reproduced only as algebra from
  (3),(4)+(10),(11); no independent check of Marshall et al. (1985) was made.
- **Fig. 16 splitting-induced notch ductility, Section 3.5, Appendix kernel**:
  conceptual; the Appendix integral equation (A4) is the dislocation
  formulation implemented (in spirit) by solver 1.

## Open questions

1. Table 3's middle column is printed "2K²/E'σ₀δ₀"; we verified it equals
   2𝒢/σ₀δ₀ = 2c₁²λ (R-curve driving force with K_tip = 0) to 3 digits on all
   rows, but the K in the header is not explicitly defined in the caption.
2. The softening-strength curve of Fig. 11 was read off at 4 points (±0.03–0.05
   reading error); our values agree within that, but a tabulation (Bao & Zok
   1992, "in preparation") would allow a tighter check.
3. Eq (13): the OCR text shows "δ_b"; the printed symbol is S_b. Confidence
   high from context (units: S_bR/τ is a length).

## Extraction items of lowest confidence

- **Eq (20) exponent**: the bracket in the printed formula carries a small
  "−2" superscript; we initially misread it as +². The −² reading is forced by
  (18)+(19) and by the sentence "the tougher the matrix, the smaller the
  bridging zone", and our numerical steady-state route confirms it (a +²
  bracket would be off by a factor [(√(1+g)+√g)/(√(1+g)−√g)]²). [resolved,
  high confidence]
- **Table 1 (resolved, no erratum)**: an earlier draft transcribed the
  atomic-bond δ₀ as 10⁻⁹ m and flagged the printed σ₀δ₀ = 1 J/m² as an
  erratum. The page image prints δ₀ = 10⁻¹⁰ m; the table is internally
  consistent. The erratum claim was withdrawn and the transcription corrected.
- **Eqs (37),(38),(40)**: subscripts of f̄₂…f̄₄ partially ambiguous in print;
  only the structure (used in 39, 40) was transcribed. Not needed numerically.
- **Fig. 11 lower curve normalization**: the inset suggests both laws share
  σ₀, δ₀; our softening computation on that basis matches the drawn curve,
  which supports the reading.

## Measured error sources (test tolerances are set from these)

- Zone solver: Dugdale identities exact to 2e-5 (N = 48, L ≤ 0.05);
  Table 2 max errors |Δc₁| = 0.0040, |Δc₂| = 0.0083 (N = 128 at large β,
  where a tip boundary layer develops); Table 3 max errors |Δc₁| = 0.0020,
  |Δ2𝒢| = 0.0008.
- Spectral solver: junction log-singularity at the trailing edge limits
  openings to ~0.3–0.5 % at N = 600–1500 (measured on the eq-45 configuration,
  where the Tada closed form is available); K_tip to ~0.05 %.
- Strength search: bisection refinement on δ_t = 1; rectilinear vs eq (47)
  agreement 0.06–0.65 % at N = 700 confirms combined discretization + search
  accuracy.
