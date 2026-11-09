# ANSWER

## Results

1. **Remote stress at first connector failure: Σ\* ≈ 133.7 MPa**  (133.73 MPa)
2. **Added crack length at that instant: L\* ≈ 4.09 mm per end**  (4.086 mm; effective crack
   half-length c\* = a + L\* = 6.486 mm)

## Method

Plane strain, so E′ = E/(1−ν²) = 141.059 GPa. Model as a bridged (Dugdale-type) crack:
the separated region is |x| < c with c = a + L, remote tension Σ, and — until the first
connector fails — a uniform closing traction t = 160 MPa on a < |x| < c. Quasistatic advance
of each outer front with intrinsic work Γ₀ gives the LEFM tip condition

  K_tip = Σ√(πc) − (2tc/√(πc))·arccos(a/c) = K₀ ≡ √(E′Γ₀) = 1.7616 MPa·√m,

which fixes the branch Σ(c) = K₀/√(πc) + (2t/π)·arccos(a/c). (Check: it starts at
Σ(c→a) = K₀/√(πa) = 20.29 MPa, i.e. initiation; small L gives the textbook
small-scale-bridging zone L = (π/8)(K_applied − K₀)²/t².)

The first connector to fail is the one adjacent to the original slit, at x = a⁺, since the
opening decreases monotonically from the slit center out to the front. Its opening is the
zone-root opening. By superposition (linearity in the loads) it is the classical Dugdale
CTOD at the zero-K stress Σ_D = (2t/π)arccos(a/c) for the same (a, c), plus the Griffith
opening produced by the extra remote stress K₀/√(πc):

  δ(a) = (8ta/πE′)·ln sec θ + 4K₀√(c² − a²)/(E′√(πc)),  θ = arccos(a/c).

First failure: δ(a) = δ_c = 9.00 μm. Both terms increase monotonically with c (from 0 at
c = a), so the root is unique; brentq gives c\* = 6.4856 mm, hence L\* = 4.0856 mm and
Σ\* = Σ(c\*) = 133.73 MPa. This is genuine large-scale bridging (L\* ≈ 1.7a), so Σ\* sits well
below the small-scale steady-state estimate √(E′(Γ₀ + tδ_c)/(πa)) = 165.4 MPa.

## Checks

- **Dugdale CTOD anchor:** integrating the singular-integral (dislocation-density) solution
  of the same elastic problem reproduces (8ta/πE′)ln sec θ to 6×10⁻¹¹ relative error, and
  reproduces δ(a) = 9.000000 μm at (Σ\*, c\*).
- **Independent solver:** a completely separate spectral (Chebyshev) solution of the
  hypersingular integral equation, solving the full 2×2 system (K_tip = K₀, δ(a) = δ_c)
  without the closed form, gives Σ = 133.70 MPa, L = 4.079 mm — within 0.02% / 0.17%
  (its truncation error, which alternates with resolution and brackets the closed-form value).
- **Profile shape:** the computed opening is maximum at the slit center and decreases
  monotonically to zero at the front (δ(0)/δ(a) ≈ 1.53), confirming the connector at x = a⁺
  fails first; K_tip and δ → 0 as c → a (initiation consistency).
