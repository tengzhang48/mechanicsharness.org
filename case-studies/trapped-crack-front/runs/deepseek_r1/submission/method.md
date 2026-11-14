# Method

## 1. Geometry, loading and the first-order front-perturbation kernel

The crack is planar (`y = 0`) with front `x = a(z)`. With

* `n = (1, 0, -a')/sqrt(1+a'^2)` (outward normal in the crack plane),
* `t = n x e_y = (a', 0, 1)/sqrt(1+a'^2)`,

the local SIFs are `K1 = lim sqrt(2 pi r) sigma_yy`,
`K2 = lim sqrt(2 pi r) e_y.sigma.n`, `K3 = lim sqrt(2 pi r) e_y.sigma.t`.
For a straight front `n = e_x`, `t = e_z`, and the prescribed loading gives
`K1^0`, `K2^0 = 0.6 K1^0`, `K3^0 = 0.8 K1^0`, uniform along the front.

For a small in-plane front perturbation `da(z)` the standard first-order
(Rice 1985) result for a planar crack in an infinite body is

```
dK_a/K_a^0 = (1/(2 pi)) Int [ da(z') - da(z) ] / (z'-z)^2 dz' ,   a=1,2,3.
```

The kernel is the same for the three modes (the near-tip scalar problems for
the relevant displacement component have the same mixed Neumann/Dirichlet
structure).  For a Fourier mode `da = A exp(i k z)` its symbol is

```
dK_a/K_a^0 = -(1/2) |k| da ,      |k| = 2 pi / lambda .
```

The sign is fixed by stability/thermodynamics: an advanced part of the front
has a smaller SIF (and a retarded part a larger one), so a straight front is
stable under mode I and the front is pushed back where the toughness is high.
This is the convention used below.  The geometric rotation of `(n,t)` is a
separate first-order effect and is added explicitly.

## 2. Part A: prescribed cosine front

Write `a(z) = abar + A cos(kz)`, `k = 2 pi/lambda`, and `eps_A = A/lambda`.
The nonlocal kernel gives the common relative perturbation

```
F(z) = -(1/2) k A cos(kz) = -pi eps_A cos(kz).
```

The rotation of the local frame (to first order in `a' = -A k sin(kz)`) gives

```
K2 = K2^0 - a' K3^0 + ... ,
K3 = K3^0 + a' K2^0 + ... ,
K1 = K1^0 (1 + F) + O(eps_A^2) .
```

Therefore

```
K1 = K1^0 + eps_A K1^0 [ -pi cos ] ,
K2 = K2^0 + eps_A K1^0 [ -0.6 pi cos + 1.6 pi sin ] ,
K3 = K3^0 + eps_A K1^0 [ -0.8 pi cos - 1.2 pi sin ] .
```

For the energy release rate, with
`G = (1-nu^2)(K1^2+K2^2)/E + (1+nu)K3^2/E` and
`G^0 = [(1-nu^2)(1+0.6^2) + (1+nu)0.8^2] (K1^0)^2/E`, one obtains

```
G/G^0 = 1 + eps_A [ -2 pi cos - 2 pi beta sin ] + O(eps_A^2),
beta = 2[ -(4/3) w2 + (3/4) w3 ],
```

where `w1,w2,w3` are the normalized energy weights of the three modes.  For
`nu = 0.3`, `beta = 0.1809045226...`.  Thus the eight Part A numbers are

```
C1 = -pi,       S1 = 0,
C2 = -0.6 pi,   S2 = 1.6 pi,
C3 = -0.8 pi,   S3 = -1.2 pi,
CG = -2 pi,     SG = -2 pi beta.
```

## 3. Part B: equilibrium front in the tough bands

Let `a(z) = abar + eps lambda f(s)`, `s = z/lambda`, and
`K_a = K_a^0 + eps K1^0 kappa_a(s)`.  The toughness is

```
Gamma(z) = Gamma0 [1 + eps g(s)],   g(s) = (1-q^2)/(1-2q cos 2 pi s+q^2) - 1
        = 2 sum_{n>=1} q^n cos(2 pi n s),      q = 1/2.
```

The nonlocal kernel gives the relative perturbation `F = -H[f']/2`, where
`H` is the Hilbert transform (`H[f']` has Fourier symbol `2 pi |n|`).  The
geometric rotation gives the same `-0.8 f'` and `+0.6 f'` terms as in Part A.
Hence

```
kappa1 = -H[f']/2 ,
kappa2 = 0.6 kappa1 - 0.8 f' ,
kappa3 = 0.8 kappa1 + 0.6 f' .
```

The energy release rate to first order is

```
G/G^0 = 1 + eps [ 2 kappa1 + beta f' ] = 1 + eps [ -H[f'] + beta f' ] .
```

Equilibrium `G = Gamma` at every point therefore reduces to

```
-H[f'] + beta f' = g(s).
```

Its Fourier solution (for `g_n = q^n`, `n >= 1`) is

```
f_n = q^n (-1 - i beta) / (2 pi n (1+beta^2)),
```

or, in closed form,

```
f(s) = [ (1/2) ln D(s) - beta phi(s) ] / (pi (1+beta^2)),
D(s) = 1 - 2 q cos(2 pi s) + q^2,
phi(s) = arg[1 - q exp(2 pi i s)] = atan2(-q sin(2 pi s), 1-q cos(2 pi s)).
```

The first term is the familiar symmetric crack-trapping logarithm; the
`beta phi` term is an antisymmetric contribution caused by the mixed-mode
(mode II + mode III) loading and the rotation of the local shear basis.  It is
why the front is not symmetric about the tough-band centre and why the
maximum of `kappa3` is not at `s = 0` or `s = 1/2`.

### Loading magnitude

Averaging `G = Gamma` over a period removes all first-order perturbation
terms (both `-H[f']` and `beta f'` have zero mean), so in the limit
`eps -> 0`

```
G^0 = Gamma0
=> (K1^0)^2 [(1-nu^2)(1+0.6^2) + (1+nu)0.8^2] / E = Gamma0
=> K1^0 sqrt((1-nu^2)/(E Gamma0)) = sqrt((1-nu^2)/2.0696) = 0.6630976493...
```

### Extrema

Maximizing/minimizing the closed-form `f` and `kappa3` over one period gives

```
max f - min f = 0.3436236624...   (min at s = -0.01429903, max at s = 0.45733145)
kappa3_max    = 1.2744082244...   at s = 0.0448009387...
```

The 16 requested samples `s = k/16`, `k = 0,...,15`, are listed in
`results.json`; they are samples of the continuous functions, not the
locations of the extrema (the peak-to-peak and `kappa3_max` above are the
continuous extrema over the period).

## 4. Checks performed

* The closed-form Fourier coefficients of `f` reproduce the numerical FFT of
  the closed form to `~2e-16`.
* The equilibrium equation `-H[f'] + beta f' = g` is satisfied to
  `~6e-13` (limited by the FFT resolution; the analytic identity is exact).
* Substituting the first-order `kappa_a` into the Irwin relation reproduces
  `G/G^0 - 1 = eps g(s)` at all sampled points to machine precision (the
  residual scales to zero with `eps`).
* The Part A coefficients were also checked by directly evaluating
  `G(z)` from the perturbed `K_a(z)` and fitting the `cos`/`sin` amplitudes;
  the result agrees with `CG` and `SG` to machine precision.
* The sign of the kernel was checked against the thermodynamic requirement
  that a localized front bump costs positive line energy (equivalently, an
  advanced part of the front has a lower SIF).  This fixes the sign of
  `C1 = -pi` and makes the front recede at the tough bands (`f(0) < 0`).

## 5. Assumptions and caveats

* First order in the front slope (`A/lambda` in Part A, `eps` in Part B);
  all terms `O(eps^2)` are dropped.
* The Rice kernel is used in its standard infinite-body, planar-crack form
  with coefficient `1/(2 pi)` and is assumed to act identically on modes
  I, II and III.  This mode universality follows from the scalar
  Neumann/Dirichlet structure of the three near-tip fields, but it is the
  main external input of the calculation.
* The local SIFs are defined with the local normal `n` and tangent `t`; the
  geometric rotation terms `-a' K3^0` and `+a' K2^0` are included.  Their
  rotation acting on the nonlocal perturbation is second order and dropped.
* The loading is the one that is uniform for every straight front; only the
  zeroth-order mean `G^0 = Gamma0` fixes its magnitude (there is no `O(eps)`
  shift of the loading because a uniform front translation does not change
  the SIFs).

## 6. Files

* `solve.py` -- computes the numbers analytically and writes `results.json`;
  includes the Fourier and equilibrium checks.
* `results.json` -- the requested deliverable.
