# A crack front trapped by tough bands under mixed-mode loading

Solve the problem below and report the requested numbers in `results.json`. The problem is fully stated here; the
method is your choice. Python 3 with NumPy and SciPy is available. Work only inside this folder.

## 1. The body and the crack

- An infinite, homogeneous, isotropic, linear elastic solid with Young's modulus E and Poisson's ratio **ν = 0.3**.
  Small strains, quasi-static loading, no body forces.
- Cartesian coordinates (x, y, z). A planar crack lies in the plane y = 0 and occupies the region x < a(z), where
  a(z) is a single-valued, smooth function. The crack front is the curve x = a(z), y = 0, which extends to
  z = ±∞. The crack faces are traction-free and do not touch.
- The crack can extend only within the plane y = 0 (it runs along a weak plane of the solid).

**Stress intensity factors.** At a point of the front, let **n** be the unit vector normal to the front, lying in
the plane y = 0 and pointing away from the crack (towards x > a), and let **t** = **n** × **e**<sub>y</sub>. For a
straight front x = const this gives **n** = **e**<sub>x</sub> and **t** = **e**<sub>z</sub>. With σ the stress
tensor evaluated in the plane y = 0 at distance r from the front along **n**, the stress intensity factors are

- K<sub>1</sub> = lim<sub>r→0</sub> √(2πr) σ<sub>yy</sub>,
- K<sub>2</sub> = lim<sub>r→0</sub> √(2πr) **e**<sub>y</sub>·σ·**n**,
- K<sub>3</sub> = lim<sub>r→0</sub> √(2πr) **e**<sub>y</sub>·σ·**t**.

**Energy release rate** per unit area of crack advance at a point of the front (Irwin's relation):

G = (1 − ν²)(K<sub>1</sub>² + K<sub>2</sub>²)/E + (1 + ν) K<sub>3</sub>²/E.

## 2. The loading

The crack is loaded from infinity only, and the loading does not change when the front moves. Its defining
property is this: for any straight front a(z) = const, wherever it lies, the stress intensity factors are uniform
along the front and equal to

K<sub>1</sub><sup>0</sup> > 0,   K<sub>2</sub><sup>0</sup> = 0.6 K<sub>1</sub><sup>0</sup>,   K<sub>3</sub><sup>0</sup> = 0.8 K<sub>1</sub><sup>0</sup>.

(Equivalently: every front considered here stays within a bounded distance of a straight line, and at distances
from the front, measured across it in a plane z = const, that are large compared with the wavelength and the
amplitude of its shape, the stresses tend to the leading, inverse-square-root stresses of a straight front with
these three factors.) G<sup>0</sup> denotes the energy release rate of a straight front.

## 3. Part A: a prescribed wavy front

Let the front be a(z) = ā + A cos(2πz/λ) with λ > 0 and A/λ ≪ 1 (ā is any constant). To first order in A/λ the
stress intensity factors and the energy release rate along this front can be written

K<sub>α</sub>(z) = K<sub>α</sub><sup>0</sup> + (A/λ) K<sub>1</sub><sup>0</sup> [C<sub>α</sub> cos(2πz/λ) + S<sub>α</sub> sin(2πz/λ)],   α = 1, 2, 3,

G(z) = G<sup>0</sup> + (A/λ) G<sup>0</sup> [C<sub>G</sub> cos(2πz/λ) + S<sub>G</sub> sin(2πz/λ)].

Report the eight dimensionless numbers C<sub>1</sub>, S<sub>1</sub>, C<sub>2</sub>, S<sub>2</sub>, C<sub>3</sub>, S<sub>3</sub>,
C<sub>G</sub>, S<sub>G</sub> for ν = 0.3 and the loading of Section 2.

## 4. Part B: the equilibrium front in tough bands

The fracture energy of the weak plane varies along z only:

Γ(z) = Γ<sub>0</sub> [1 + ε g(z/λ)],   g(s) = (1 − q²)/(1 − 2q cos 2πs + q²) − 1,   **q = 0.5**,

with 0 < ε ≪ 1. The function g has period 1 in s and zero mean; the tough bands are centred at s = 0, ±1, ±2, ….

The front is in equilibrium: **G(z) = Γ(z) at every point of the front**. The magnitude of the loading (at the
fixed ratios of Section 2) is the one for which such an equilibrium front exists.

Work to first order in ε. Write the front and the stress intensity factors as

a(z) = ā + ε λ f(z/λ),   with f bounded, periodic with period 1 (f(s + 1) = f(s)) and of zero mean over a period,

K<sub>α</sub>(z) = K<sub>α</sub><sup>0</sup> + ε K<sub>1</sub><sup>0</sup> κ<sub>α</sub>(z/λ),   α = 1, 2, 3,

where K<sub>α</sub><sup>0</sup> are the factors of the straight front under the same loading. Report:

1. `K10_scaled` = K<sub>1</sub><sup>0</sup> √((1 − ν²)/(E Γ<sub>0</sub>)) in the limit ε → 0;
2. `front_peak_to_peak` = max f − min f;
3. f(s), κ<sub>1</sub>(s), κ<sub>2</sub>(s) and κ<sub>3</sub>(s) at the sixteen points s = k/16, k = 0, 1, …, 15;
4. `kappa3_max`, the largest value of κ<sub>3</sub>(s) over a period, and `s_kappa3_max`, the point in
   −1/2 ≤ s < 1/2 where it occurs.

## 5. What to deliver

- `results.json` in this folder, in exactly this form (numbers as JSON numbers; lists in the order k = 0, …, 15):

```
{
  "partA": {"C1": 0.0, "S1": 0.0, "C2": 0.0, "S2": 0.0, "C3": 0.0, "S3": 0.0, "CG": 0.0, "SG": 0.0},
  "partB": {"K10_scaled": 0.0, "front_peak_to_peak": 0.0,
            "f": [16 numbers], "kappa1": [16 numbers], "kappa2": [16 numbers], "kappa3": [16 numbers],
            "kappa3_max": 0.0, "s_kappa3_max": 0.0}
}
```

- `method.md`: how you obtained the numbers, the assumptions you made, and the checks you ran.
- Any code you used (including for checks), in this folder.

Each reported number is compared with a reference value; it passes if it is within 1% of the reference value's
magnitude or within 0.005, whichever is larger.
