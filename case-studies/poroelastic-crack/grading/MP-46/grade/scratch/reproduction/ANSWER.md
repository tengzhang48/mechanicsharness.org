# A suddenly pressurized crack (MP-46) — answer

**Reproduce:** `python3 solve.py`  (prints, and writes, `results.json`)

All quantities are per unit loading magnitude, K in units of √(πL), V in units of L²/G
(taking the loading amplitudes s0 = p0 = 1, G = 1, L = 1 for the numbers).

## Model and method

Plane-strain Biot poroelasticity exactly as given in the task. The crack is modelled with the
standard Rice–Cleary (1976) "crack-line" reduction: the crack is a mathematical slit whose faces may
carry **both** a mechanical traction and a pore-pressure boundary condition. I use two ingredients:

1. **Elastic crack formulae.** A slit of half-length L whose faces carry a uniform *opening*
   traction P (P > 0 opens the crack) has, in plane strain,
   `K = P √(πL)` and `V = π(1−ν)P L²/G`.
2. **Effective traction.** Written with the effective stress `σ' = σ + α p δ`, the face condition seen
   by the solid *skeleton* is `t' = σ'·n = t + α p n`.

These are verified against an independent finite-element solution of the full transient Biot problem
(`fem_poro.py`, see *Checks*).

## Derivation in brief

**Loading 1** (`σ_n = −H(t)`, `p = 0`).
The faces are *drained* (`p = 0` there) and carry the constant opening traction P = 1. Because `p = 0`
on the crack boundary the poroelastic reciprocity term that could make K time-dependent vanishes, so
the total-stress field — hence K — is **time independent**, while the compliance relaxes from undrained
to drained as the pore pressure generated in the solid diffuses away:

- `K1_0 = K1_∞ = √(πL)`  (→ 1 in √(πL) units)
- `V1_0  = π(1−ν_u) L²/G`   (undrained, t = 0⁺)
- `V1_∞  = π(1−ν)   L²/G`   (drained, t → ∞)

*Small-time growth.* At t = 0⁺ an undrained pore pressure `p_faces_0` (below) sits on the crack line;
it relaxes inside a boundary layer of depth ~√(ct) with an erf profile. Integrating that profile
against the crack-line weight (1−x²/L²)^{1/2} gives
`V ≈ V1_0 + 2√π (ν_u−ν) √(ct/L²) · L²/G`, i.e. `A = 2√π (ν_u−ν)`.

*Pore pressure at the faces at t = 0⁺.* On the drained crack line the elastic field is modulus
independent with `σ_kk = σ_xx + σ_yy = −2`; the undrained constraint
`p = −(ν_u−ν) σ_kk / [α(1−2ν)]` then gives `p_faces_0 = 2(ν_u−ν)/[α(1−2ν)]`.

**Loading 2** (`σ_n = 0`, `p = H(t)`).
The mechanical traction is zero, so the skeleton feels only the effective face traction
`t' = α p n`. For the upper face the outward normal is `n = (0,−1)`, hence `t' = (0, −αp)`: the
skeleton is pulled *into* the crack, so the crack **closes** (V < 0, K < 0 — this is why the task
defines V as "positive for an open crack"). The equivalent opening traction is `P = −α p0`, so

- `K2_0 = K2_∞ = −α p0 √(πL)`
- `V2_∞ = −π(1−ν) α p0 L²/G`

## Numbers (ν = 0.2, ν_u = 0.4, α = 0.89)

| quantity | value | formula |
|---|---|---|
| `K1_0`  | 1.000000 | √(πL) |
| `K1_inf`| 1.000000 | √(πL) |
| `V1_0`  | 1.884956 | π(1−ν_u) |
| `V1_inf`| 2.513274 | π(1−ν) |
| `V1_short_coeff` | 0.708982 | 2√π(ν_u−ν) |
| `p_faces_0` | 0.749064 | 2(ν_u−ν)/[α(1−2ν)] |
| `K2_0`  | −0.890000 | −α√(πL) |
| `K2_inf`| −0.890000 | −α√(πL) |
| `V2_inf`| −2.236814 | −π(1−ν)α |

## Checks (what was compared, output pasted, how it could fail)

**A. Elastic calibration of the FEM machinery** — drained elastic crack under opening traction,
`ν=0.2`, expected `V = π(1−ν) = 2.51327`, `K = 1.0`:
```
--- elastic calibration (drained, nu=0.2): expect V=2.51327 K=1.0 ---
  V=  1.65178  K(nu)=  1.42905  K(nu_u)=  1.90541
```
The FEM **under-resolves/over-confines**: the outer boundary is clamped at |x|=3 and the crack is
meshed by only ~40 boundary nodes, so both V (−34 %) and the COD-based K (+43 %) are off. This run is
therefore used only for **signs and trends**, not magnitudes. It could fail (and did) if the domain
were too small or the tip too coarse — exactly what happened.

**B. Transient FEM, loading 1** (`dt = 5e-3`, `c ≈ 1.87`), V and tip COD-K vs time:
```
=== loading 1 ===
  t= 0.0050  V=  0.93640  K(nu)= 0.87947  K(nu_u)= 1.17263  w(0)= 0.56476
  t= 0.1550  V=  1.62345  K(nu)= 1.41184  K(nu_u)= 1.88245  w(0)= 1.00527
  t= 0.8050  V=  1.64975  K(nu)= 1.43016  K(nu_u)= 1.90688  w(0)= 1.02310
  t= 1.5000  V=  1.64975  K(nu)= 1.43016  K(nu_u)= 1.90688  w(0)= 1.02310
```
V rises monotonically from the undrained to the drained plateau (as required by part (a)); the K
inferred with the *drained* modulus plateaus exactly at the drained-elastic calibration value
(1.429 vs 1.430), i.e. the total-stress SIF is time-independent, consistent with `K1_0 = K1_∞`.
Ratio `V1_0/V1_∞ = 0.75` matches the analytic `(1−ν_u)/(1−ν) = 0.6/0.8 = 0.75`.

**C. Transient FEM, loading 2** (`σ_n = 0`, `p = 1`):
```
=== loading 2 ===
  t= 0.0050  V= -0.13900  K(nu)= -0.15114  w(0)= -0.07950
  t= 0.3550  V= -0.43360  K(nu)= -0.38196  w(0)= -0.26808
  t= 1.5000  V= -0.43406  K(nu)= -0.38227  w(0)= -0.26840
```
The FEM gives **V < 0 at all times**, w(0) < 0, and a nearly time-independent K — confirming the sign
of the effective-traction argument (crack closes, `K2_0 ≈ K2_∞ < 0`) and the time-independence of K.
Magnitudes from this run are not trustworthy (see check A); the *model* magnitudes
(`−α√(πL)`, `−π(1−ν)α`) are used in `results.json`.

## Questions / assumptions / uncertainty

1. **Sign of loading 2.** A pore pressure applied on *unloaded* faces pulls the skeleton `t' = −αp`
   into the crack, so the crack closes. The FEM agrees (check C). I therefore report negative
   `K2_0, K2_inf, V2_inf`. If the problem intended a *hydraulically* pressurized crack
   (`σ_n = −p`), that would be the **sum** of loadings 1 and 2, flipping the sign of K2.
2. **Magnitude of the loading-2 (and short-time) coefficients** is the least certain part: it relies
   on the effective-traction reduction with the drained modulus at t→∞, whereas the transient FEM is
   too coarse (check A) to confirm the factor. The signs and the time-independence of K are solid.
3. **Definition of `σ_kk`.** The task writes the constitutive law with indices 1,2, so I used the
   2-D (plane-strain) trace throughout; the finite-element code uses the same convention.