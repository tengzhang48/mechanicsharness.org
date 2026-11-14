# Method

## 1. Coplanar crack-front perturbation

I used the first-order perturbation for a small **coplanar** advance
`delta a(z)` of a straight crack front.  In Fourier space, with front-shape
wavenumber `k` and `beta = 2(1-nu)/(2-nu)`,

```
dK1_hat = -|k| K10 delta_a_hat / 2
dK2_hat = -|k| K20 delta_a_hat / 2 + i k K30 delta_a_hat
dK3_hat = -|k| K30 delta_a_hat / 2 + beta i k K20 delta_a_hat.
```

Multiplication by `i k` is `d/dz` in physical space.  The choice of `t` in
the problem fixes the signs of the shear-coupling terms.  For the given
loading, `K20 = 0.6 K10` and `K30 = 0.8 K10`.

For Part A, `delta a = A cos(kz)` and `k lambda = 2 pi`.  Substitution gives

```
C1 = -pi,                 S1 = 0
C2 = -0.6 pi,            S2 = 1.6 pi
C3 = -0.8 pi,            S3 = -1.2 beta pi.
```

Define

```
D = (1-nu^2)(1 + 0.6^2) + (1+nu) 0.8^2 = 2.0696.
```

Thus `G0 = K10^2 D/E`.  Directly differentiating Irwin's relation gives
`CG = -2 pi` and `SG = -0.46803508643489666`.  In particular, the sine term
in `G` comes from the mixed-mode II/III coupling; it is not obtained by
assuming that each SIF has only a cosine perturbation.

## 2. Loading and equilibrium front

At zeroth order in `epsilon`, equilibrium and the zero mean of `g` require
`G0 = Gamma0`.  Therefore

```
K10_scaled = sqrt((1-nu^2)/D) = 0.6630976492661631.
```

Let `H = |D_s|`, the periodic operator with Fourier symbol `2 pi |n|`, and
let a prime denote `d/ds`.  The energy perturbation can be written

```
delta G/G0 = epsilon (-H f + eta f') + O(epsilon^2),
eta = 2(0.6)(0.8) [beta(1+nu) - (1-nu^2)] / D
    = 0.0744900975465563.
```

Consequently the first-order equilibrium equation is

```
H f - eta f' = -g,       mean(f) = 0.
```

The toughness profile has the Fourier coefficients

```
g_n = q^|n|  (n != 0),    g_0 = 0,
g(s) = 2 sum_{n>=1} q^n cos(2 pi n s).
```

Solving mode by mode gives

```
f_n = -q^|n| / [2 pi (|n| - i eta n)].
```

Equivalently, the real series used for evaluation is

```
f(s) = -1/[pi(1+eta^2)]
       sum_{n>=1} q^n/n [cos(2 pi n s) - eta sin(2 pi n s)].
```

The SIF functions follow without another equilibrium solve:

```
kappa1 = -H f/2
kappa2 = -0.6 H f/2 - 0.8 f'
kappa3 = -0.8 H f/2 + 0.6 beta f'.
```

I evaluated the complex Fourier series with `|n| <= 96`; since `q=0.5`,
the omitted tail is negligible at floating-point precision.  The sixteen
samples are exactly those at `s=k/16`.  Global extrema were bracketed on a
dense periodic grid and refined with `scipy.optimize.minimize_scalar`.
The refined front minimum and maximum are respectively

```
f_min = -0.219856662487379 at s = 0.994079103672757
 f_max =  0.128789923558966 at s = 0.482253671936630,
```

which gives the reported peak-to-peak value.  The same procedure was used
for `kappa3_max` and its location.

## 3. Code and checks

- `solve.py` derives all coefficients, solves the Fourier system, constructs
  `results.json`, refines the extrema, and prints numerical residuals.
- `solve_checks.txt` records the solver checks.  On 20,000 points per
  period, the maximum equilibrium residual, direct energy residual,
  reconstructed-energy residual, and Fourier-versus-closed-form `g`
  residual were all `2.22e-15` or smaller.
- `verify.py` independently evaluates the solution using real closed-form
  sums (logarithms and arctangents), rather than the solver's complex
  series.  It checks every JSON sample, the equilibrium identity
  `H f - eta f' = -g`, the energy identity reconstructed from all three
  `kappa` arrays, the loading normalization, `CG = -2 pi`,
  `SG = -2 pi eta`, and the refined extrema.  All these checks pass.

Run with:

```
python3 solve.py
python3 verify.py
```

The calculations are consistent to first order in both `A/lambda` and
`epsilon`.  Higher-order terms and a separate numerical solution of the
full three-dimensional elasticity problem were not computed.
