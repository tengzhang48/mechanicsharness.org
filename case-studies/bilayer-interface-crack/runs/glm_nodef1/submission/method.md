# Method

## Problem setup and normalization

Units: `h = 1`, `H = 2` (h/H = 0.5), `E2 = 1`, `E1 = 9`, ν1 = ν2 = 1/3, plane strain;
P = 1 (case A), M = 1 (case B). Reported quantities are dimensionless and
unaffected by this choice. Plane-strain moduli Ēi = Ei/(1−ν²) = 10.125 and 1.125.
Dundurs parameters (material 1 = upper layer): α = 0.8, β = 0.2, and the
oscillation index ε = ln[(1−β)/(1+β)]/(2π) = ln(2/3)/(2π) = −0.06453178.
The crack-tip field is written in the standard form

    T(r) = σ22 + i σ12 = K r^{iε} / sqrt(2πr)   (r ahead of the tip),

with **r measured in units of h**, so A = |K| and, at r = h (r/h = 1), arg T = arg K.
All phase angles below are arg K under this convention (see "Assumptions").

## Energy release rate G (exact)

Far ahead of the tip (x1 → +∞) the bonded bilayer is stress-free (the layer-2
stresses given in the task have zero resultant force and moment together with
those of layer 1). Far behind the tip each layer is in its exact Saint-Venant
state (the given σ11 distributions satisfy equilibrium and the traction-free
outer faces exactly, so they are exact elastic solutions, not approximations).
The J-integral evaluated far from the tip therefore gives exactly the difference
of far-field strain-energy densities:

    G = Σi [ Ni²/(2 Ēi ti) + 6 Mi²/(Ēi ti³) ]

with arm resultants (N1, M1) = (−P, −M), (N2, M2) = (P, M + P(h+H)/2):
- case A: G E2 h/P² = 1/16.2·… = **1.7716049383** (terms 0.0493827 + 0.2222222 + 1.5)
- case B: G E2 h³/M² = **1.2592592593** (terms 0.5925926 + 0.6666667)

These are the Suo–Hutchinson "beam" values and are exact, not beam approximations.

## G–K relation (derived and verified)

From a Riemann–Hilbert analysis of the crack-tip eigenpair (potentials
F1 = a z^λ, F2 = −a z^λ, Ω1′ = −ā c* z^λ̄, Ω2′ = +ā c* z^λ̄, λ = −1/2+iε_m,
c* = e^{−2πiλ̄}; `code/eigenpair.py` verifies traction continuity ahead,
displacement continuity ahead, traction-free crack faces, and reproduces the
FEM-measured traction and crack-opening amplitudes), the J-integral of the
eigenpair integrated numerically gives, to 5×10⁻⁷ relative,

    G = |K|² (1/Ē1 + 1/Ē2) / (2 cosh²(πε))  =  |K|² · 64/135 · E2⁻¹,

i.e. sech²(πε) = 24/25 exactly for this material pair. The homogeneous limit
reproduces G = |K|²/Ē exactly. Hence
- case A: A√h/P = |K_A| = sqrt(G_A/c) = **1.9331268**
- case B: A h^{3/2}/M = |K_B| = **1.6298006**

## Phases and ε (numerical, Q8 FEM)

`code/fem.py`: 8-node quadrilateral plane-strain FEM on a structured mesh,
graded geometrically toward the crack tip (tip elements down to 4×10⁻⁶ h in x1
and 10⁻⁵ h in x2), outer surfaces traction-free, end x1 = +L traction-free,
end x1 = −L (L = 30h) loaded with the exact Saint-Venant stress profiles as
consistent nodal forces; interface nodes duplicated for x1 < 0. Loads are
self-equilibrated, so no artificial boundary conditions are needed (3 rigid
modes pinned at a far corner; L-independence verified: L = 15/30/45 change |K|
by <0.02% and the phase by <0.001°).

Extraction (`code/fits.py`): K from a least-squares fit of
T(r)√(2πr) r^{−iε} = K + c√r + d r over several windows (r/h ∈ [5×10⁻⁴, 4×10⁻²]),
and independently from the crack-opening fit δ(r) = Q K̄ r^{1/2−iε}(1 + …) with
the eigenpair-derived Q. Six mesh-refinement levels per load case
(`code/final_run.py`, audit in `code/final_audit.txt`):

- phase A: values across the finest meshes/windows/both fits span 52.25–52.30°;
  reported (consensus median) **52.27°**
- phase B: span 4.24–4.26°; reported **4.25°**
- |K| from FEM converges monotonically to the exact-relation values above
  (finest mesh: 1.9301 and 1.6273, i.e. −0.16% and −0.15% — consistent within
  the observed convergence trend).
- ε itself was measured from the FEM log-slopes of arg[T√(2πr)] (≈ −0.0645)
  and of arg[δ/√r] (≈ +0.0645, opposite sign), confirming the analytic
  ε = −0.06453178 (magnitude and sign).
- phase change per decade = ε ln 10 = **−8.5136°**.
- combined loading M = Ph: superposition K_tot = K_A + K_B was verified against
  a direct (P = M = 1) FEM run (agreement to 6 digits); with the exact |K|s and
  the phases above, arg K_tot = **30.43°**.

## Checks run (all passed)

1. Element self-tests: Q8 shape-function derivatives vs finite differences
   (1×10⁻¹⁰), patch test (uniform strain to 1×10⁻¹⁰), K·1 = 0.
2. Homogeneous benchmarks with exact answers, same pipeline:
   - moments ±M on equal arms: exact K_I = 2√3 M h^{−3/2}; FEM 3.4634
     (−0.02%), phase −0.014° (exact 0) — also validates the √(2π)
     normalization and the phase extraction end-to-end.
   - forces ∓P on equal arms: exact |K| = √7 P h^{−1/2} (the arm-2 moment
     P(h+H)/2 = P makes this ≠ P); FEM 2.6460 (+0.01%).
3. G–K relation vs closed form S/(2cosh²πε): 5×10⁻⁷ (two material pairs).
4. T-fit vs δ-fit K: agree to 0.1% in magnitude, 0.05° in phase.
5. FEM |K| vs exact-relation |K|: agree to 0.2% (slow monotone convergence).
6. Superposition check; domain-length independence; solver residual < 10⁻⁴.
7. Displacement-jump amplitude matches the eigenpair constant Q (0.1%).

## Assumptions

- The phase at r = h is arg[K (r/h)^{iε}] at r/h = 1, i.e. **arg K** — the
  only unit-independent choice (h is the problem's length scale). A different
  convention (e.g. absorbing (2πr)^{iε} into the power) would shift all three
  reported phase angles by the same constant (e.g. by ε ln 2π ≈ −6.8°).
- "Ignore contact" is taken literally: the oscillatory near-tip interpenetration
  is accepted; quantities are defined by the crack-tip field continued from r → 0.
- Loads applied at x1 = −30h (Saint-Venant decay of end detail makes this
  equivalent to x1 → −∞; verified by the L study).
- The sign convention for T(r), P (compressive in layer 1), and M follows the
  stress formulas in the task statement.

## What remains unverified

- The r = h phase convention (assumption above); all three phase numbers shift
  together under an alternative convention.
- An additional end-to-end bimaterial benchmark (central interface crack under
  remote tension, Rice–Sih configuration) was attempted but abandoned: my
  custom two-tip mesh generator produced numerically singular systems (cause
  not fully diagnosed). It is not used anywhere in the reported quantities.
- The FEM phases carry a possible systematic bias of order 0.01–0.05°
  (estimated from the homogeneous benchmark, where the exact phase is 0).

## Files

- `results.json` — the deliverable.
- `code/fem.py` — Q8 FEM (mesh, assembly, loads, solves, evaluation).
- `code/fits.py` — K/ε extraction and fitting utilities.
- `code/eigenpair.py`, `code/check_c.py` — crack-tip eigenpair, J-integral,
  G–K relation.
- `code/final_run.py` — production run producing `results.json`
  (audit: `code/final_audit.txt`).
- `code/selftest.py`, `code/bench_homog.py`, `code/validations.py`,
  `code/stab1.py`, `code/stab2.py`, `code/run_bimaterial.py` — checks.
- `code/debug*.py`, `code/central_crack.py`, `code/central2.py` — debugging
  trail of the abandoned benchmark (kept for transparency).
