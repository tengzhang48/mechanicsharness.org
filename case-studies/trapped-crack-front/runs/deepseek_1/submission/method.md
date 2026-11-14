# Method

## Summary

I used the first-order 3-D crack-front perturbation theory of Rice (1985), in the form in
which the local stress-intensity-factor (SIF) change produced by a small in-plane front
perturbation `δa(z)` is written as

```
δK_α(z) = K_β^0 [ -A_{αβ} H[δa](z) + B_{αβ} (dδa/dz) ],
```

where `H` is the non-local operator whose Fourier symbol is `|k|` (`H e^{ikz} = |k| e^{ikz}`).
The first term is the "advance" (front-position) term; the second is the local "tilt" term.
The unperturbed straight front has `K_1^0, K_2^0 = 0.6 K_1^0, K_3^0 = 0.8 K_1^0`.

### Response matrices used

In the mode basis (1,2,3) = (I,II,III):

* `A = diag(1/2, -1/2, -1/2)`. That is, a forward bump lowers `K_I` (mode I is stable) and
  raises `K_II` and `K_III` (modes II and III are unstable). This is the standard Rice
  result; the II/III degeneracy is the shear-vector degeneracy of the isotropic problem.
* `B` has only `B_23 = -1` and `B_32 = +1`. These follow directly from the geometric
  projection of the reference SIF vector onto the local frame `(n, e_y, t)` for a front
  `x = a(z)`: `n = (1,0,-a')/√(1+a'^2)`, `t = n × e_y = (a',0,1)/√(1+a'^2)`, hence
  `K_II^loc = K_II - a' K_III` and `K_III^loc = K_III + a' K_II`. There is no I–II tilt
  coupling for a smooth in-plane front: the crack plane remains `y = 0`, so the mode-I
  direction `e_y` is not rotated.

The G response follows from Irwin's relation
`G = (1-ν²)(K_1²+K_2²)/E + (1+ν)K_3²/E`, linearised about the straight front.

## Part A

For `a(z) = ā + A cos(2πz/λ)`, let `k = 2π/λ`, so `a' = -A k sin(kz)`. Substituting into
the response law gives

* `C_α = -2π A_{αα} (K_α^0/K_1^0)`,
* `S_α = -2π B_{αβ} (K_β^0/K_1^0)`,

and for the energy release rate

* `C_G = -2π C_adv`, `S_G = -2π C_tilt`,

with

```
C_adv  = 2 [ w (A_11 K_1² + A_22 K_2²) + v A_33 K_3² ] / G^0
C_tilt = 2 K_2 K_3 ( w B_23 + v B_32 ) / G^0
w = (1-ν²)/E,   v = (1+ν)/E
G^0 = w (K_1² + K_2²) + v K_3² .
```

Numerically (ν = 0.3, K_2/K_1 = 0.6, K_3/K_1 = 0.8, E = 1):

```
G^0/K_1² = 2.0696
C_adv = -0.1206030...,   C_tilt = +0.1809045...
```

## Part B

### Loading magnitude

At ε = 0 the equilibrium condition is `G^0 = Γ_0`. With the above `G^0`,

```
K_1^0 √((1-ν²)/(E Γ_0)) = √((1-ν²) / 2.0696) = 0.66309765...
```

### Equilibrium front

Write `a(z) = ā + ε λ f(s)`, `s = z/λ`. The toughness is
`Γ = Γ_0[1 + ε g(s)]` with the Poisson-kernel function

```
g(s) = (1-q²)/(1-2q cos 2πs + q²) - 1 = 2 Σ_{n≥1} q^n cos(2π n s),   q = 1/2.
```

To first order the equilibrium condition `G = Γ` is

```
g(s) = -C_adv H[f](s) + C_tilt f'(s).
```

In Fourier space (`f = Σ [f_n^c cos(2πns) + f_n^s sin(2πns)]`):

```
f_n^c = -C_adv q^n / (π n Δ²)
f_n^s =  C_tilt q^n / (π n Δ²)
Δ² = C_adv² + C_tilt²
```

whose closed form is

```
f(s) = (1/(πΔ²)) [ (C_adv/2) ln(1 - 2q cos 2πs + q²)
                   + C_tilt arctan( q sin 2πs / (1 - q cos 2πs) ) ].
```

This is the mixed-mode generalisation of the log-sine front shape found by Gao & Rice for
mode-I crack trapping. The `C_tilt` term makes the front asymmetric. Its derivative, used
below, is

```
f'(s) = [2q/(Δ² D)] [ C_adv sin(2πs) + C_tilt (cos(2πs) - q) ],
D = 1 - 2q cos(2πs) + q² .
```

The SIF perturbations are then obtained from

```
κ_α(s) = (K_α^0/K_1^0) [ -A_{αα} H[f](s) ] + Σ_β (K_β^0/K_1^0) B_{αβ} f'(s),
```

with

```
H[f] = (C_tilt f' - g)/C_adv .
```

The reported `f`, `κ_1`, `κ_2`, `κ_3` are these functions at `s = k/16`, `k = 0,...,15`.
`front_peak_to_peak` is `max f - min f` over one period; `kappa3_max`/`s_kappa3_max` are the
maximum of `κ_3` and its location in `[-1/2, 1/2)`.

## Checks performed

1. **Energy consistency.** At all sixteen sample points, the reported `κ_α` satisfy
   `2[w(K_1 κ_1 + K_2 κ_2) + v K_3 κ_3]/G^0 = g(s)` to machine precision. This is a
   nontrivial check of the algebra connecting the SIF response to the G response.
2. **Fourier vs closed form.** The front shape and all `κ_α` were evaluated both by direct
   Fourier summation (4000–8000 harmonics) and from the closed form; agreement was
   `~1e-15`.
3. **Mode-I limit.** Setting `K_2^0 = K_3^0 = 0` reduces the front equation to
   `g = -H[f]` with `C_adv = 1`, giving
   `f = (1/2π) ln(1 - 2q cos 2πs + q²)`, the standard Gao–Rice log-sine trapping shape.
4. **2-D limits of the crack operator.** The Fourier symbols used (mode I/II:
   `μ|k|/[2(1-ν)]`, mode III: `μ|k|/2`) reduce correctly to the classical 2-D
   plane-strain and antiplane crack operators.
5. **G balance in Part B.** The reported `κ_3` turns out to be exactly `3.316... g(s)`
   for this loading (the 3-4-5 mix makes the `f'` contributions cancel), so `κ_3` peaks at
   the tough-band centre `s = 0`.

## Assumptions and caveats

* The calculation is first order in `ε` and in `A/λ`; no higher-order or finite-amplitude
  corrections are included.
* The Rice (1985) response matrices are used. The mode-I coefficient `-1/2` and the
  II–III geometric tilt are well established; the equal, opposite-sign treatment of modes
  II and III (`+1/2` for both) is the standard "mode I stable, shear modes unstable"
  result but is the main unverified assumption of this solution. A different mode-II
  coefficient would change `C_2`, `C_G`, `f`, `κ_1` and `κ_2` (but not `K10_scaled`).
* The loading is taken to be the self-similar K-field loading of Section 2; no additional
  uniform-stress terms are included.
* The front is assumed smooth and single-valued; the tilt term is the local linear
  projection, and no sharp-kink (Cotterell–Rice) I–II coupling is included because the
  crack plane stays `y = 0`.
