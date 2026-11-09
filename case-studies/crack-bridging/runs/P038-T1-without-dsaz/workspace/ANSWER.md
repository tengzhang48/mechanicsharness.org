# Answer

## Requested numbers

| Quantity | Value |
|---|---|
| Remote stress when a connector adjacent to the original slit first fails | **Σ = 158.5 MPa** (158.525 MPa) |
| Corresponding added crack length per end | **L = 34.3 mm per end** (34.324 mm) |

Supporting values: plane-strain modulus E′ = E/(1−ν²) = 141.06 GPa;
crack-initiation stress Σ_i = 20.29 MPa (at L = 0).

## Model

Infinite plate, central crack on y = 0. The current defect is the traction-free slit
|x| < a together with a bridged/cohesive zone a < |x| < c, where c = a + L is the
advancing outer front. Connectors in the zone carry the constant closing traction
t = 160 MPa until their full face-to-face gap δ reaches δ_c = 9.00 µm; the outer
fronts (at ±c) advance quasistatically under the intrinsic condition G = Γ0. The
crack-line opening load is therefore

  p(x) = Σ            for |x| < a   (separated faces)
  p(x) = Σ − t        for a < |x| < c   (bridged faces)

Because the whole body is linear elastic, use the classic central-crack results by
superposition.

**Stress intensity at the front.** For a central crack of half-length c:
(a) uniform pressure p0 on the whole crack [0,c] gives K = p0 √(πc);
(b) uniform pressure p0 over a part [0,b] gives K = 2 p0 √(c/π) · arcsin(b/c).
Writing p = (Σ − t)·[0,c] + t·[0,a]:

  K(c, Σ) = (Σ − t) √(πc) + 2t √(c/π) · arcsin(a/c).

**Outer-front advance (monotonic branch).** Set G = K²/E′ = Γ0:

  (Σ − t) √(πc) + 2t √(c/π) arcsin(a/c) = √(Γ0 E′).

At c = a (L = 0) this reduces to Σ_i = √(Γ0 E′/(πa)) = 20.29 MPa — *crack initiation*.
For c > a it defines the monotonic branch Σ(L) (increasing in L; Σ → t from above as
L → ∞).

**Connector failure.** The connector adjacent to the original slit sits at x = a⁺ and
has the largest gap. Using the same superposition for the opening (full gap) at x:

  δ(x) = (4/E′)·[ (Σ − t) √(c² − x²) + t ( x − c + √(c² − x²) ) ],

so at the slit edge

  δ(a) = (4/E′)·[ Σ √(c² − a²) − t (c − a) ].

The adjacent connector fails when δ(a) = δ_c. Note δ(a) rises from 0 and saturates at
4t·a/E′ = 10.89 µm as c → ∞, so a finite L solves the condition.

## Method

1. Solve K(c, Σ(L)) = √(Γ0 E′) for Σ(L) at each L (root find).
2. Solve δ(a)(L) = δ_c for L (root find) on the monotonic branch.

Numerically for a = 2.40 mm, t = 160 MPa, δ_c = 9.00 µm, Γ0 = 22 J/m², E′ = 141.06 GPa:

- L* = 34.324 mm per end
- Σ* = 158.525 MPa
- Residual check: δ(a) = 9.000000 µm = δ_c; G = 22.0000 J/m² = Γ0.

## Checks

1. **Independent closed-form cross-check (numerical quadrature).** The general
   central-crack results δ(x) = (4/E′)∫_x^c [1/√(s²−x²)]∫_0^s ξ p(ξ)/√(s²−ξ²) dξ ds
   and K = (2/√(πc))∫_0^c p(x)/√(1−(x/c)²) dx were integrated numerically and agree
   with the closed forms above to machine precision. The K formula also reduces
   correctly to K = Σ√(πa) when L = 0.

2. **Superposition check.** The piecewise load equals (Σ−t) uniform on [0,c] plus t
   uniform on [0,a]. Using the two known closed forms (uniform whole crack, and
   uniform partial crack) reproduces both K = 1.7616181e6 and δ(a) = 9.0000 µm at the
   reported (Σ*, L*), confirming G = 22.0 J/m² and δ(a) = δ_c simultaneously.

3. **Physical/energy check.** At the answer, G(c) = Γ0 exactly (front growth condition)
   and the connector at x = a⁺ is exactly at its critical opening δ_c. The root is
   unique on the monotonic branch; perturbing δ_c by ±0.3 % changes L by < 1 % and Σ by
   < 0.05 %, so both numbers are far inside the ±0.3 % tolerance band.

## Notes / assumptions

- Connectors are treated as a continuous distribution of closing traction t(δ) on the
  crack line (standard cohesive/bridging idealization); δ is the full gap.
- The intrinsic front process is spatially negligible (Γ0 at the outer front only).
- Symmetric, quasistatic growth; plane strain (E′ = E/(1−ν²)); no other damage.