# Answer: mode III (anti-plane) power-law crack-tip field, amplitude fixed by J

## Requested numbers

| # | Quantity | Value |
|---|---|---|
| 1 | Δw(r = 0.0200 mm) | **1.318 × 10⁻³ mm** (≈ 1.32 µm, upper face +z) |
| 1 | Δw(r = 0.0800 mm) | **1.794 × 10⁻³ mm** (≈ 1.79 µm) |
| 2 | τ_rz(r = 0.0200 mm, θ = 60°) | **+38.59 MPa** |
| 2 | τ_θz(r = 0.0200 mm, θ = 60°) | **+138.01 MPa** |
| 3 | r where \|γ\| = 0.0100 on θ = 0° | **0.0469 mm** |
| 3 | r where \|γ\| = 0.0100 on θ = 90° | **0.0251 mm** |
| 3 | r where \|γ\| = 0.0100 on θ = 180° (upper face, limit) | **0.0134 mm** |

(Signed per the positive e_r, e_θ basis; the \|γ\| = 0.01 level corresponds to τ = τ₀(0.01/γ₀)^{1/n} = 126.71 MPa. Check: Δw ratio 1.794/1.318 = 1.361 = (0.08/0.02)^{1/(n+1)} ✓.)

## Model chosen and why

The problem is **anti-plane (mode III)**: only w(x,y) ≠ 0, w odd in y, faces at θ = ±π traction-free (τ_θz = 0 there). The given law γ = (γ₀/τ₀)(τ/τ₀)^{n−1}τ is a pure power law, i.e. a nonlinear-elastic (deformation-theory) material; since loading is monotonic with a fixed pattern, stresses at a fixed point in the dominated zone stay proportional, so deformation theory is exact here, a strain-energy density exists, and J is path-independent (knowledge/TOPIC.md §2, §6.6).

The leading field is the separable mode III power-law (Hult–McClintock/HRR-type) solution, stress function ψ = a r^m g(θ):

- τ_rz = a r^{m−1} g′(θ),  τ_θz = −a m r^{m−1} g(θ),  ρ²(θ) = g′² + m² g²;
- **stress exponent fixed by J**: a separable field has J ∝ r^{1−p(n+1)}, so a finite nonzero J = 0.120 N/mm forces p = 1/(n+1) = 0.2222, i.e. m = n/(n+1) = 0.7778, τ ~ r^{−1/4.5}, γ ~ r^{−3.5/4.5}, w ~ r^{1/4.5} (TOPIC.md §3);
- angular function g from compatibility (curl-free γ), a 2nd-order ODE (ρ^{n−1}g′)′ + mα ρ^{n−1} g = 0, α = 1/(n+1), with BCs g′(0) = 0 (antisymmetry, τ_rz(0) = 0), g(π) = 0 (traction-free faces), normalization g(0) = −1 (solved by shooting, `solve.py`);
- displacement w = (γ₀ aⁿ /(α τ₀ⁿ)) r^α ρ^{n−1}(θ) g′(θ), with no rigid translation (w → 0 as r → 0).

**How the loading fixes the magnitude:** the separable field has a single free amplitude a, and J = C_J a^{n+1} with C_J = 2.8466 × 10⁻¹⁰ (computed from the angular integral with the g(0) = −1 normalization). Setting J = 0.120 N/mm gives **a = 82.534 MPa·mm^p**; equivalently τ ∝ J^{1/(n+1)}, Δw ∝ J^{n/(n+1)}. All requested points (r ≤ 0.08 mm) lie well inside the r < 0.30 mm zone where the leading field is stated to dominate.

## Assumptions added/changed

None of the problem's assumptions were changed. I adopted (i) the deformation-theory reading of the given total stress–strain law (exact under the stated monotonic, fixed-pattern loading) and (ii) the leading-term-only field, as invited by the problem ("use that leading field"). Everything else (τ₀, γ₀, n, J, traction-free faces, antisymmetry, no elastic term) is taken exactly as given.

## Verification (`solve.py`, all pass)

1. **Eigen-consistency:** with m = n/(n+1) imposed by J-independence, the face condition residual g(π) = 2 × 10⁻¹⁷; independently shooting m to enforce g(π) = 0 recovers m = 0.777777777778 = n/(n+1) to 12 digits.
2. **Energy flux (main quantitative check):** J recomputed from the scaled field by brute-force quadrature on circular arcs at r = 0.020, 0.080, 0.150 mm, with w,ₓ obtained by Cartesian finite differences of w(x,y) (avoiding the polar-basis error flagged in TOPIC.md §4): **0.120000 N/mm at all three radii**.
3. **Known limiting solution (n = 1):** the solver reproduces g = −cos(θ/2) to 2 × 10⁻¹⁵ and the amplitude relation reduces to J = K²/(2μ), μ = τ₀/γ₀ = 40 000 MPa, with a = K√(2/π), to 1.3 × 10⁻¹⁵ relative.
4. Field residuals: equilibrium 7.6 × 10⁻⁷, compatibility 1.5 × 10⁻⁷ (relative, finite differences); τ ∥ γ to 10⁻¹⁶; face traction τ_θz(π) = −3 × 10⁻¹⁵ MPa; w > 0 on the upper face as stated.

## Missing / uncertain

- Higher-order (non-singular) terms from the outer geometry are neglected per the problem (corrections called negligible below 0.30 mm); their size is not quantifiable from the given data.
- Small strain fails at the very tip (γ → ∞ as r → 0); |γ| = 1 only inside r ≈ 2 × 10⁻⁴ mm, far below all requested radii, so the numbers above are unaffected.
- The pure power law has no elastic term, so stresses diverge (τ ~ r^{−1/4.5}); this is the specified law, but real materials would cut this off.
- Results scale as J^{1/(n+1)} (stresses) and J^{n/(n+1)} (Δw), so any error in the given J propagates fractionally by those factors.

## Do these results determine the load at slit advance?

**No.** The field and its J-value are loading-response (amplitude) quantities; onset of growth requires a separate material fracture criterion (e.g. a critical J or critical face slip Δw), which is not contained in the given data or in these results.
