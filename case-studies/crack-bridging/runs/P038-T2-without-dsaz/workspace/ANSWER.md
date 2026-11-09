# Separation between unequal strips — Answer

## Results

| Quantity | Value |
|---|---|
| 1. Applied moment per strip | **M = 0.2126 N·m** |
| 2. Front-to-root distance | **L = 10.80 mm** |
| 3. Opening at mid-zone | **δ(L/2) = 38.22 µm** |

(Exact computed values: M = 0.21262579 N·m, L = 10.8001 mm, δ = 38.2189 µm.)

---

## Model and justification

**Geometry / loading.** Two strips of width b = 25 mm bonded along a plane. Upper:
E₁ = 110 GPa, h₁ = 0.60 mm. Lower: E₂ = 190 GPa, h₂ = 0.90 mm. A 20 mm starter
cut is opened at the left ends. Equal and opposite *opening couples* M (one per
strip, no end forces) are applied at the free ends and increased quasistatically.

**Flexural rigidities (full width).**

    EI₁ = E₁ b h₁³/12 = 110e9 · 0.025 · (0.6e-3)³/12 = 4.950e-2 N·m²
    EI₂ = E₂ b h₂³/12 = 190e9 · 0.025 · (0.9e-3)³/12 = 2.8856e-1 N·m²
    R  = 1/EI₁ + 1/EI₂ = 23.6675 (N·m²)⁻¹

Axial stretch and shear are neglected; each arm is an Euler–Bernoulli beam and the
full relative opening is δ = w₁ + w₂.

**Why pure moments make this a steady-state (front-independent) problem.** With a
couple and no end force, the bending moment in each separated arm equals the
applied M (no transverse load outside the cohesive zone). For two unequal arms the
relative-opening curvature in the open part is δ″ = M R, and the available energy
release rate is

    G = M² R / (2b),

which is **independent of crack length** (the hallmark of the pure-moment DCB).

**Cohesive (retained-connector) zone.** Behind/around the advancing front the new
faces retain connectors that pull the faces shut. Since the closing traction is a
bolt load on each arm, the relative opening obeys the beam-on-elastic-foundation
equation

    δ⁗ = − b·t(δ)·R ,    t(δ) = σ_b (1 − δ/δ_c).

With σ_b = 0.20 MPa, δ_c = 0.20 mm this is

    δ⁗ − k δ = − b σ_b R ,   k = b σ_b R/δ_c = 5.9169e8 m⁻⁴ ,  β = k^{1/4} = 155.96 m⁻¹

(oscillatory/decaying, wavelength 2π/β ≈ 40.3 mm). The particular solution is
δ_p = δ_c, so

    δ(x) = δ_c + A cosh βx + B sinh βx + C cos βx + D sin βx.

**Boundary conditions.**

* At the **front** (x = 0) the crack tip has just failed: the connector there is at
  δ = δ_c (traction zero), and it matches the open crack behind it, where the arm
  moment is M and the shear is zero:

      δ(0) = δ_c ,   δ″(0) = M R ,   δ‴(0) = 0 .

* At the **root** (x = L, the far end of the retained/process zone) the faces come
  back together and join the long straight bonded region:

      δ(L) = 0 ,   δ′(L) = 0 .

The event "a newly retained connector first reaches δ_c" is exactly the tip
condition δ(0) = δ_c: the front has advanced until the trailing connector of the
zone is at the failure opening, i.e. the cohesive/retained zone is fully formed
and spans from the front back to the original starter-cut root.

**Moment from the energy (R-curve) condition.** Advancing the front costs
Γ₀ = 1.40 J/m² of intrinsic fracture energy *plus* the connector work. The area
under the linear-softening law is ½ σ_b δ_c = ½·0.20e6·0.20e-3 = 20.0 J/m²
(per unit joint area, pair of faces counted once). Hence the crack-growth
resistance — and therefore the applied G — is

    G = Γ₀ + ½ σ_b δ_c = 1.40 + 20.0 = 21.40 J/m².

Solving G = M² R/(2b):

    M = √( 2 b G / R ) = √( 2·0.025·21.40 / 23.6675 ) = 0.21263 N·m.

**Length L and mid-opening.** Using this M, the two conditions at x = L give
(eliminating B) a single scalar equation whose smallest physically admissible root
is L = 10.80 mm. The opening at the midpoint is then

    δ(L/2) = δ_c + A cosh(βL/2) + B sinh(βL/2) + C cos(βL/2) + D sin(βL/2)
           = 38.22 µm .

Because of the oscillatory kernel, δ decays from δ_c = 200 µm at the front to 0
at the root over ~10.8 mm; at the middle the value is 38.2 µm.

---

## Checks

1. **All boundary conditions satisfied:** computed δ(0) = 200.000 µm (= δ_c),
   δ″(0) = 5.0323 = M R, δ‴(0) = 0, δ(L) = 0 (−1.1e-9 residual), δ′(L) = 0
   (−3.7e-18). The 4th-order ODE residual at L/2 is 0 (analytic), confirming the
   profile is an exact solution of δ⁗ − kδ = − b σ_b R.

2. **Energy-release formula cross-check:** the compliance method for the two-arm
   pure-moment member gives dθ/da = M R, so G = ½ M (dθ/da)/b = M²R/(2b) = 4.7335e2
   · M² — identical to the closed form used. With M obtained from G = Γ₀ + ½σ_bδ_c,
   the supplied energy release equals the total resistance (stable, front-independent
   growth criterion).

3. **Isolated root:** the scalar length equation has an isolated zero at
   L = 10.800 mm (residual 1e-9) between sign changes at 8 mm (+5.8e-5) and 12 mm
   (−6.9e-6), confirming a unique process-zone length and ruling out spurious roots.

---

## Scripts

`solve.py` (main model), `verify.py` (BC + ODE + energy checks),
`check.py` (G-formula verification). All use only NumPy/SciPy.