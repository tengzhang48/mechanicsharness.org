# An interface crack between two elastic layers

Solve the problem below and report the requested numbers in `results.json`. The problem is fully stated here; the
method is your choice. Python 3 with NumPy and SciPy is available. Work only inside this folder.

## 1. The bilayer and the crack

- Plane strain. Two isotropic, linear elastic layers, each infinitely long in the x<sub>1</sub> direction, are bonded
  along the interface x<sub>2</sub> = 0. Small strains, no body forces.
- **Layer 1** (above, 0 < x<sub>2</sub> < h): Young's modulus E<sub>1</sub> = 9 E<sub>2</sub>, Poisson's ratio
  ν<sub>1</sub> = 1/3.
- **Layer 2** (below, −H < x<sub>2</sub> < 0): Young's modulus E<sub>2</sub>, Poisson's ratio ν<sub>2</sub> = 1/3.
- Thickness ratio **h/H = 0.5**.
- The outer surfaces x<sub>2</sub> = h and x<sub>2</sub> = −H are traction-free.
- A semi-infinite crack lies on the interface along x<sub>1</sub> < 0, with its tip at the origin. The layers are
  perfectly bonded for x<sub>1</sub> > 0. The crack faces are traction-free; ignore any contact between them.

## 2. The loading

The bilayer is loaded only far behind the crack tip (x<sub>1</sub> → −∞), where each layer is a separate beam. There
the stresses are σ<sub>12</sub> = σ<sub>22</sub> = 0 and

- in layer 1:   σ<sub>11</sub> = −P/h − 12 M (x<sub>2</sub> − h/2)/h³,
- in layer 2:   σ<sub>11</sub> = P/H + 12 (M + P (h + H)/2)(x<sub>2</sub> + H/2)/H³.

So layer 1 carries an axial compressive force P and a bending moment M (positive M bends layer 1 away from layer 2),
and the stresses in layer 2 are those that leave the bonded part far ahead of the tip (x<sub>1</sub> → +∞)
stress-free. Forces and moments are per unit length in the x<sub>3</sub> direction.

- **Load case A:** P > 0, M = 0.
- **Load case B:** P = 0, M > 0.

## 3. The crack-tip field

Write the tractions on the interface a distance r ahead of the tip (x<sub>1</sub> = r > 0, x<sub>2</sub> = 0) as one
complex number, T(r) = σ<sub>22</sub> + i σ<sub>12</sub>. As r → 0 the tractions are dominated by the singular
crack-tip field. The quantities asked below belong to that crack-tip field. Where a value at a finite distance is
asked, it means the crack-tip field continued to that distance, not the full stress there.

G is the energy release rate (energy released per unit area of crack advance along the interface).

## 4. What to report

1. Load case A: `G_E2_h_over_P2` = G E<sub>2</sub> h/P²; `amplitude_sqrt_h_over_P` = A √h/P, where
   A = lim<sub>r→0</sub> √(2πr) |T(r)|; `phase_at_h_deg`: the phase angle arg T of the crack-tip field at r = h.
2. Load case B: `G_E2_h3_over_M2` = G E<sub>2</sub> h³/M²; `amplitude_h32_over_M` = A h<sup>3/2</sup>/M;
   `phase_at_h_deg`, as in case A.
3. `phase_at_h_both_deg`: the same phase angle when both loads act together with M = P h (P > 0).
4. `phase_change_per_decade_deg`: for load case A, the change of the phase angle arg T of the crack-tip field when r
   increases tenfold. Report 0 if it does not change.

Phase angles in degrees, −180° < angle ≤ 180°.

## 5. What to deliver

- `results.json` in this folder, in exactly this form (numbers as JSON numbers):

```
{
  "caseA": {"G_E2_h_over_P2": 0.0, "amplitude_sqrt_h_over_P": 0.0, "phase_at_h_deg": 0.0},
  "caseB": {"G_E2_h3_over_M2": 0.0, "amplitude_h32_over_M": 0.0, "phase_at_h_deg": 0.0},
  "phase_at_h_both_deg": 0.0,
  "phase_change_per_decade_deg": 0.0
}
```

- `method.md`: how you obtained the numbers, the assumptions you made, and the checks you ran.
- Any code you used (including for checks), in this folder.

Each reported number is compared with a reference value. The energy release rates and the amplitudes pass if they are
within 1% of the reference value; each phase angle at r = h passes if it is within 0.5° of the reference value; the
change per decade passes if it is within 0.2° of the reference value.
