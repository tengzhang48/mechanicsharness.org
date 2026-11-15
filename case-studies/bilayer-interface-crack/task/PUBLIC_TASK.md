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

## 3. Definitions

With μ<sub>k</sub> = E<sub>k</sub>/(2(1 + ν<sub>k</sub>)) and κ<sub>k</sub> = 3 − 4ν<sub>k</sub> (plane strain), the
bimaterial constant is

ε = (1/2π) ln[(κ<sub>1</sub>/μ<sub>1</sub> + 1/μ<sub>2</sub>)/(κ<sub>2</sub>/μ<sub>2</sub> + 1/μ<sub>1</sub>)].

The complex stress intensity factor K = K<sub>1</sub> + i K<sub>2</sub> is defined by the tractions on the interface a
distance r ahead of the tip (x<sub>1</sub> = r > 0, x<sub>2</sub> = 0):

σ<sub>22</sub> + i σ<sub>12</sub> = K r<sup>iε</sup> / √(2πr)   as r → 0,

where r<sup>iε</sup> = exp(iε ln r). The phase angle of the loading at the reference length h is

ψ = arg(K h<sup>iε</sup>),   in degrees, −180° < ψ ≤ 180°.

(ψ does not depend on the unit of length, as long as r and h are measured in the same unit.)

G is the energy release rate (energy released per unit area of crack advance along the interface).

## 4. What to report

1. `eps`: the bimaterial constant ε.
2. Load case A: `G_E2_h_over_P2` = G E<sub>2</sub> h/P²; `Kabs_sqrt_h_over_P` = |K| √h/P; `psi_deg` = ψ.
3. Load case B: `G_E2_h3_over_M2` = G E<sub>2</sub> h³/M²; `Kabs_h32_over_M` = |K| h<sup>3/2</sup>/M; `psi_deg` = ψ.
4. `psi_both_deg`: the phase angle ψ when both loads act together with M = P h (P > 0).

## 5. What to deliver

- `results.json` in this folder, in exactly this form (numbers as JSON numbers):

```
{
  "eps": 0.0,
  "caseA": {"G_E2_h_over_P2": 0.0, "Kabs_sqrt_h_over_P": 0.0, "psi_deg": 0.0},
  "caseB": {"G_E2_h3_over_M2": 0.0, "Kabs_h32_over_M": 0.0, "psi_deg": 0.0},
  "psi_both_deg": 0.0
}
```

- `method.md`: how you obtained the numbers, the assumptions you made, and the checks you ran.
- Any code you used (including for checks), in this folder.

Each reported number is compared with a reference value. `eps`, the energy release rates and the magnitudes of K pass
if they are within 1% of the reference value; each phase angle (`psi_deg`, `psi_both_deg`) passes if it is within 0.5°
of the reference value.
