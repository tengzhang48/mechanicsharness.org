# Problem: a suddenly pressurized crack in a fluid-saturated porous solid

An infinite, homogeneous, isotropic, fluid-saturated porous solid deforms in plane strain
(x–y plane). A crack occupies −L ≤ x ≤ L, y = 0. The solid obeys Biot's linear poroelasticity, with
shear modulus G, drained Poisson's ratio ν, undrained Poisson's ratio ν_u (ν < ν_u ≤ ½), Biot
coefficient α (0 < α ≤ 1) and permeability coefficient κ.

In terms of the total stress σ_ij (tension positive), the pore pressure p, the strain e_ij and the
variation of fluid content per unit volume ζ:

- 2G e_ij = σ_ij − ν σ_kk δ_ij + α(1 − 2ν) δ_ij p
- 2G ζ = α(1 − 2ν) σ_kk + α²(1 − 2ν)² p / (ν_u − ν)
- equilibrium: σ_ij,j = 0
- Darcy: q_i = −κ p,i
- mass balance: ∂ζ/∂t + q_i,i = 0

(indices 1, 2; summation implied). The diffusivity is c = 2κG(1 − ν)(ν_u − ν) / [α²(1 − 2ν)²(1 − ν_u)].

Stresses and pore pressure are zero everywhere for t < 0 and vanish at infinity. Two loadings act on
both crack faces (−L ≤ x ≤ L) for t > 0. Here σ_n = σ_yy is the normal stress on the faces and H is
the unit step.

- **Loading 1:** σ_n = −H(t), p = 0.
- **Loading 2:** σ_n = 0, p = H(t).

Let K be the mode I stress intensity factor at the crack tips, and V the crack volume per unit
thickness (the integral of the opening over the crack; positive for an open crack).

(a) For loading 1, find K and V just after loading (t = 0⁺) and at long times (t → ∞). Find also the
leading-order growth of V at small times, ct/L² ≪ 1.

(b) For loading 2, find K at t = 0⁺ and as t → ∞, and V as t → ∞.

(c) Evaluate your results for ν = 0.2, ν_u = 0.4 and α = 0.89.

**Report in `results.json`** (for (c); K in units of √(πL), V in units of L²/G, per unit loading
magnitude):
- `K1_0`, `K1_inf`, `V1_0`, `V1_inf`: loading 1 at t = 0⁺ and t → ∞;
- `V1_short_coeff`: the coefficient A in V ≈ V1_0 + A √(ct/L²) for loading 1 at small times;
- `p_faces_0`: in loading 1, the pore pressure induced in the solid immediately next to the crack
  faces at t = 0⁺, before it diffuses;
- `K2_0`, `K2_inf`, `V2_inf`: loading 2;
- `formulas`: your general expressions for all of the above, as strings.
