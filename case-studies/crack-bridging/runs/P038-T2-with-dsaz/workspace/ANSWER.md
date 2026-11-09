# Separation between unequal strips — answer

## Requested numbers (state: first connector reaches δc and goes inactive)

| quantity | value |
|---|---|
| opening couple `M` on each strip | **0.2126 N·m** |
| bridging-zone length `L` | **10.80 mm** |
| opening at mid-zone `δ(L/2)` | **38.22 µm** |

(Full precision: M = 0.2126258 N·m, L = 10.80010 mm, δ(L/2) = 38.21890 µm.)

## Mechanical description

The specimen is a **double-cantilever-beam (DCB) splitting specimen** with two
unequal arms bonded far ahead:

* free ends (left) → **starter cut** `a0 = 20 mm` (no connectors) → **root**
  (`x = L`) → **bridging zone** (`0 < x < L`, connectors) → **front**
  (`x = 0`) → long, straight, **bonded** region far ahead;
* `x` runs from the front (0) toward the root (L); free ends at `x = L + a0`.

Each strip is an Euler–Bernoulli beam (stretching and shear neglected as
stated). The two strips are **clamped at the front**: the bonded region far
ahead is long and stays straight, so it neither rotates nor translates —
exactly the standard built-in-end condition for a DCB arm. The relative opening
is the sum of both arm deflections, `δ = w₁ + w₂`, giving the combined
rigidity `D̃ = D₁D₂/(D₁+D₂)` with `Dᵢ = Eᵢ b hᵢ³/12`:

* `D₁ = 0.049500 N·m²` (upper, E=110 GPa, h=0.60 mm)
* `D₂ = 0.288563 N·m²` (lower, E=190 GPa, h=0.90 mm)
* `D̃  = 0.042252 N·m²`

Governing equation on the bridging zone, with the given softening law
`t(δ) = σb(1 − δ/δc)`, `σb = 0.20 MPa`, `δc = 0.20 mm`:

```
D̃ δ'''' = − b·t(δ)          (0 < x < L)      connectors (closing)
D̃ δ'''' = 0                 (L < x < L+a0)   starter cut (free)
δ(0) = δ'(0) = 0                              clamped at the front (bonded part)
D̃ δ''(L+a0) = M , δ'''(L+a0) = 0          pure opening couple M on each arm
```

The connector tractions act in the bridging zone only; the starter cut is
traction free and merely transmits the constant moment `M` to the free ends.

### Fracture / front condition

With a pure opening couple and no end forces, the global (far-contour)
J-integral of the DCB is **length independent**:

```
G_far = M² / (2 b D̃) .
```

Path independence of J across the bridging zone (Rice; eq. 15 of the source
unit) gives the tip value

```
G_tip = G_far − ∫₀^{δt} t(δ) dδ ,     δt = δ(L),
```

and quasistatic advance of the front requires `G_tip = Γ₀`.  Since δ(x)
increases monotonically from the front to the root, the connector that reaches
δc first is the trailing one at the root, so the state of interest is

```
δt = δ(L) = δc .
```

At that instant every connector is still within (0, δc), so
`∫₀^{δc} t dδ = σb δc / 2`, and

```
M² / (2 b D̃) = Γ₀ + σb δc / 2 = 1.40 + 20.0 = 21.40 J/m²
```

so the moment is fixed independent of the zone length (pure-moment DCB):

```
M = sqrt( 2 · 0.025 · 0.042252 · 21.40 ) = 0.212626 N·m .
```

This is exactly the **steady-state (plateau) moment** of the R-curve
(G_ss = Γ₀ + σbδc/2): the moment plateaus precisely as the trailing connector
reaches δc.  `L` then follows from `δ(L) = δc` using the beam solution above.

## Method

1. Solve the linear ODE `D̃ δ'''' = −b t(δ)` on `0<x<L` and `δ''''=0` on
   `L<x<L+a0` analytically.  The particular solution of the softening law is
   the constant `δc`; the homogeneous solutions are `e^{±kx}, cos kx, sin kx`
   with `k = (bσb/D̃δc)^{1/4} = 1/Lc`, `Lc = 6.4118 mm`.  Eight unknowns
   (4 per region) are fixed by the clamped, continuity and free-end conditions
   (`solve.py`).
2. Trace the whole R-curve by continuation in `L`: for each `L` solve
   `G_tip = Γ₀` for `M`, then report `δ(L)`.  `M` rises smoothly and plateaus
   at `0.212626 N·m`; `δ(L)` rises monotonically and crosses `δc = 0.20 mm` at
   `L = 10.800 mm` (`trace_rcurve.py`, `final_state.py`).
3. Read `δ(L/2)` from the same solution: **38.2189 µm**.

## Checks

* **Independent shooting integration** of the ODE with `scipy.integrate`
  (`verify.py`) reproduces the analytic profile to 5 digits at `(L*, M*)`:
  δ(L)=200.000 µm, δ(L/2)=38.2189 µm, δ(free end)=2074.14 µm.
* **Independent closed form** (complex-exponential construction, `check2.py`)
  gives M = 0.21262579 N·m, L = 10.800102 mm, δ(L/2) = 38.21890 µm — identical.
* **Energy / J balance** (`verify.py`): G_far = 21.400 J/m² and
  G_tip = G_far − ∫₀^{δc}t dδ = 21.400 − 20.000 = 1.400 J/m² = Γ₀. ✓
* **Consistency with the knowledge unit**: the result is the softening-triangle
  plateau of Bao & Suo (the trailing edge reaching δ₀ marks steady state), and
  the moment is length-independent as required for pure-moment DCB loading.
* **Robustness**: `L*` and `δ(L/2)` are unchanged (to 6 digits) when the
  starter cut `a0` is varied from 5 mm to 100 mm — as they must be, since the
  moment is fixed by `G_far` and the bridging-zone solution depends only on the
  front condition.

## Files

* `solve.py` — analytic beam DCB cohesive-zone solver (main result).
* `trace_rcurve.py` — R-curve continuation trace.
* `final_state.py` — precise location of the δ(L)=δc state.
* `verify.py` — independent shooting integration + J-integral energy check.
* `check2.py` — independent complex-exponential construction.