# Answer — near-tip asymptote of a tensile power-law notch

## Requested numbers

| quantity | value |
|---|---|
| **1. singularity exponent** `p` in `σθθ(r,0) ∝ r⁻ᵖ` | **p = 0.15618**  (0.156175260286) |
| **2. stress ratio** `R = σrr(r,0)/σθθ(r,0)` | **R = 0.3733**  (0.373262313312) |

(Both are read at the same small radius on the symmetry line θ = 0, opening
amplitude positive.)

Equivalently, with the Airy stress function written `Φ = K rˢ f(θ)`, the solved
eigenvalue is

    s = 2 − p = 1.843824739714,   f(0)=1,  f′(0)=0,  f″(0) = −1.263078996343,  f‴(0)=0.

## Model and boundary conditions

The tip is modelled as a sharp **V-notch** in an isotropic, pressure-insensitive,
power-law (Norton–Hoff / total-deformation) solid in **plane strain** under
**monotonic proportional** loading — exactly the constitutive law given in the
task. Near the tip the power-law (creep) strain term dominates, and the leading
symmetric (mode-I) field is sought in the separated form used by Hutchinson
(1968) for a crack:

    Φ = K rˢ f(θ),   σ_rr = Φ,r/r + Φ,θθ/r²,   σθθ = Φ,rr,   σrθ = −∂r(Φ,θ/r),
    σ_rθ = K rˢ⁻² (1−s) f′,  σθθ = K rˢ⁻² s(s−1) f,  σrr = K rˢ⁻² (s f + f″).

The material occupies **−2π/3 < θ < 2π/3**, so the notch half-angle is
α = 2π/3. The remote y-tension is symmetric about θ = 0 and drives the
symmetric opening mode, giving the standard conditions:

- **symmetry about θ = 0:** f′(0) = f‴(0) = 0 (normalisation f(0) = 1);
- **traction-free notch faces at θ = ±2π/3:** σrθ(α)=0 ⇒ f′(α)=0 and
  σθθ(α)=0 ⇒ f(α)=0.

Substituting the separated field into plane-strain equilibrium + strain
compatibility + the power-law constitutive law yields Hutchinson's fourth-order
nonlinear angular ODE for f(θ) (the same plane-strain equation implemented in the
knowledge unit's `hrr.py`). Here the exponent **s is not fixed** (as it is for a
crack); it is the **eigenvalue** selected by the two notch-face boundary
conditions. The two unknowns (s, f″(0)) are determined by requiring
f(α)=0 **and** f′(α)=0 simultaneously. The load magnitude K cancels from both
requested numbers, so only the eigenvalue problem matters.

Physical inflection of the answer: for a **crack** (α→π, zero opening) the same
ODE returns exactly s = (2n+1)/(n+1) = 1.8 for n = 4, i.e. p = 0.2. Opening the
notch to 120° lowers the stress concentration, so the singularity is **weaker**:
p drops from 0.2 to 0.1562 (s rises from 1.8 to 1.8438). The value also respects
the finite-energy constraint p < 2/(n+1) = 0.4, as required by the task's
"finite energy in any finite neighborhood".

## Method

1. Took Hutchinson's plane-strain angular ODE (knowledge unit, eq. 31).
2. Recovered the separated-field stresses from Φ by the Airy identities.
3. Shot on (s, f″(0)) with f(0)=1, f′(0)=f‴(0)=0 to satisfy f(α)=f′(α)=0 at
   α = 2π/3.
4. Read the two requested quantities:
   - p = 2 − s,
   - σθθ(0) = s(s−1)f(0) = s(s−1), σrr(0) = s f(0) + f″(0) = s + f″(0),
     so R = (s + f″(0)) / [s(s−1)].

## Verification (quantitative)

- **Crack limit reproduced exactly.** With α = π the identical solver returns
  s = (2n+1)/(n+1) to 1e-14 for n = 3, 4, 5, 9 (e.g. n=4 ⇒ 1.8000000000 vs
  1.8). This validates both the ODE and the shooting machinery.
- **Two independent integrators agree.** The eigenvalue from a fixed-step RK4
  shooter and from `scipy.integrate.solve_ivp` agree: Δs = 1.2e-14,
  Δf″(0) = 7.6e-15.
- **Step-size independence.** s changes by < 1e-12 between 1000 and 8000 steps.
- **Independent 1-D eigenvalue route** (for each s, bracket f(α)=0 with
  brentq, then root-find f′(α)=0 on the vectorized integrator) returns
  s = 1.8438247424, p = 0.1561752576 — matching the 2-D shoot to ~5e-9.
- **Plane-strain vs plane-stress distinguishability.** The same solver on the
  plane-stress angular ODE (Hutchinson eq. 13) gives p = 0.16531 for this notch,
  different from the plane-strain p = 0.15618, confirming the plane-strain form
  required by the task is the one being used (the plane-stress crack limit also
  correctly returns s = 1.8).
- **Boundary conditions verified numerically** at θ = α = 2π/3: f(α) = 7e-13 and
  f′(α) = 4e-13, so σθθ(α) and σrθ(α) vanish (traction-free faces), and
  σ̃e(0) = 0.844 ≠ 0 (nonzero opening amplitude).
- **Finite-energy bound** p = 0.1562 < 2/(n+1) = 0.4 is satisfied.

## Files

- `notch.py` — ODE, RK4 and `solve_ivp` integrators, vectorized variant.
- `eig.py` — eigenvalue continuation from the crack (α=π) to the notch (α=2π/3).
- `answer.py` — final (p, R) with the RK4/solve_ivp cross-check.
- `bccheck.py` — traction-free boundary-condition and stress-value verification.
- `vcrack.py` — crack-limit check of the ODE for n = 3,4,5,9.
- `planestress.py` — plane-stress variant (distinguishability check).
- `vnotch1d_fast.py` — independent 1-D eigenvalue route.

## Reproduce

```
python3 answer.py         # prints p and R
python3 vcrack.py         # crack-limit sanity (s -> (2n+1)/(n+1))
python3 bccheck.py        # boundary conditions at theta = 2pi/3
```