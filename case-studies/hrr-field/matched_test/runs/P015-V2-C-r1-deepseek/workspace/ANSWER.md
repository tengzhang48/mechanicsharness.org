# Answer

## 1. Requested numbers

| # | Quantity | Value |
|---|----------|-------|
| 1 | `T` | **4816 MPa** (4816.47) |
| 2a | `C_rr(60°)` | **3.942e-3 mm^1/2** |
| 2b | `C_θθ(60°)` | **6.332e-4 mm^1/2** |
| 2c | `C_γrθ(60°) = lim √r (2ε_rθ)` | **5.731e-3 mm^1/2** |
| 3 | `C_y(90°) = lim u_y(r,90°)/√r` | **9.139e-3 mm^1/2** |
| 4 | radii on θ=0 passing all three screens | **[0.0100, 0.0613] mm** (non-empty) |

The interval endpoints: lower bound `r = 0.01000 mm` (strain screen, equality), upper bound
`r = 0.06133 mm` (σ_e = 10σ_y screen, equality); `r ≤ 0.200 mm` is not binding.

## 2. Model chosen and why

**Reduction of the constitutive law.** With σ_zz = σ_xz = σ_yz = 0 this is plane stress. The last
term is classical J2 *deformation theory* with a **bilinear** uniaxial curve. In uniaxial tension
(σ_xx = σ, others 0 ⇒ σ_e = σ, s_xx = 2σ/3):

```
ε_xx = σ/E                                              for σ ≤ σ_y
ε_xx = σ/E + (1/T − 1/E)(σ − σ_y) = σ/T − (1/T − 1/E)σ_y for σ > σ_y
```

so `E` is the elastic modulus and **`T` is the post-yield tangent modulus** (dσ/dε = T). No
power-law exponent exists, so the HRR family is not applicable except in its asymptotically linear
limit.

**Near-tip state.** Ahead of the tip σ_e ~ r^(−1/2) → ∞, so `⟨1 − σ_y/σ_e⟩ → 1`: the local response
becomes *linearly elastic* with an unchanged bulk modulus `K = E/[3(1−2ν)] = 58333.3 MPa` and a
modified shear compliance

```
1/G_eff = 2[(1+ν)/E + (3/2)(1/T − 1/E)] = 3/T + (2ν−1)/E .
```

An isotropic linear material's Mode-I stress field is universal (independent of the moduli), so the
leading field is the standard `r^(−1/2)` field; only the strain↔stress map changes. Defining

```
E_eff = 9 K G_eff/(3K + G_eff),   ν_eff = (3K − 2G_eff)/(2(3K + G_eff)),
```

algebra gives the exact identity **`E_eff = T`** (verified in `compute.py`), and ν_eff = 0.486239.

**Using J and C_zz.** On θ=0, σ_xx = σ_yy = A/√r with A = K_I/√(2π), and ε_zz = −(ν_eff/E_eff)(σ_xx+σ_yy).
Hence

```
C_zz(0) = −2 (ν_eff/E_eff) K_I / √(2π)
```

J is path independent for deformation theory; on a near-tip contour the field is the effective-linear
singular field, giving plane-stress `J = K_I²/E_eff`. Eliminating K_I:

```
E_eff / ν_eff² = 2 J / (π C_zz²) = 20371.8 .
```

Solving for G_eff (with K fixed) gives G_eff = 1620.36 MPa, E_eff = T = 4816.47 MPa,
ν_eff = 0.486239. Then K_I = √(J E_eff) = 62.074 MPa·mm^1/2 and
A = K_I/√(2π) = 24.7639 MPa·mm^1/2.

**θ = 60° coefficients.** Using the polar Mode-I stress amplitudes
`S_rr = A cos(θ/2)(1+sin²(θ/2))`, `S_θθ = A cos(θ/2)(1−sin²(θ/2))`,
`S_rθ = A sin(θ/2)cos²(θ/2)` and the effective compliance:

```
C_rr   = (S_rr − ν_eff S_θθ)/E_eff = 3.9420e-3
C_θθ   = (S_θθ − ν_eff S_rr)/E_eff = 6.3318e-4
C_γrθ  = 2(1+ν_eff) S_rθ / E_eff   = 5.7311e-3
```

**C_y(90°).** From the Mode-I displacement field
`u_y = K_I/(2G_eff) √(r/2π) sin(θ/2)(κ−cosθ)`, `κ = (3−ν_eff)/(1+ν_eff)`:

```
C_y(90°) = K_I κ /(4 G_eff √π) = 0.0025√2 (3−ν_eff)/(2ν_eff) = 9.1390e-3 .
```

**Screening interval (item 4).** On θ=0, σ_e = σ_xx = A/√r and the principal strains are
`(1−ν_eff)A/(E_eff√r)` (twice) and `−2ν_eff A/(E_eff√r)`.

- `σ_e ≥ 10σ_y` ⇒ `√r ≤ A/100` ⇒ `r ≤ 0.061325 mm`.
- `max|principal strain| ≤ 0.0500` ⇒ `0.005/√r ≤ 0.05` (since 2ν_eff A/E_eff = 0.005 exactly)
  ⇒ `r ≥ 0.010000 mm`.
- `r ≤ 0.200 mm` is not binding.

Interval = **[0.01000, 0.06133] mm**.

## 3. Assumptions changed and why

- **Deformation-theory / no unloading:** the total-strain law is interpreted as a J2 deformation
  theory evaluated under the stated monotonic proportional loading; near the tip the plastic term is
  fully active. This is exactly the regime the problem specifies.
- **Plane stress with free ε_zz:** ε_zz is *not* set to zero; it is obtained from σ_zz = 0 as
  `ε_zz = −(ν_eff/E_eff)(σ_xx+σ_yy)`. This is the quantity measured by C_zz.
- **Asymptotic effective linearity:** the `σ_y` offset in the yield term is subdominant near the tip
  (it produces only bounded strains); it does not affect the `r^(−1/2)` leading field. So the leading
  amplitudes use effective (tangent) moduli, not the elastic ones.
- **`E_eff = T`:** derived, not assumed (see `compute.py`).

## 4. Verification

1. **Kinematics vs. constitutive (independent check).** The Mode-I displacement field with the
   effective moduli was differentiated numerically in polar coordinates to obtain
   `ε_rr, ε_θθ, 2ε_rθ`; the resulting coefficients are
   `3.942048e-3, 6.33178e-4, 5.731131e-3`, matching the stress→strain values to <1e-6 relative.
2. **Numerical J integral (path independence).** Recycling the reconstructed leading field, the
   contour integral `J = ∫(W N_x − t_i u_{i,x}) ds` on full circles of radius 0.001, 0.01 and 0.1 mm
   returns **0.800000 N/mm** at every radius, confirming both path independence and `J = K_I²/E_eff`.
3. **Algebraic identity.** `E_eff − T = 0` exactly (symbolic elimination, checked numerically).
4. Substituting back gives `C_zz = −0.005000`, recovering the supplied datum.

## 5. What the screens establish; crack initiation

The three screens are *defined validity checks for the leading field*, not error estimates: σ_e ≥ 10σ_y
puts the points well inside the fully yielded zone where the asymptotic field dominates; max principal
strain ≤ 0.05 keeps infinitesimal kinematics credible; r ≤ 0.2 mm is the stated range where matching
to the surrounding body is negligible. Together they delimit the annulus [0.01, 0.0613] mm where the
result is meaningful. They say nothing about fracture.

The supplied data (`E, ν, σ_y, T, J`) characterize one monotonic loading state and do **not** determine
crack initiation: the model has no damage, cohesive, or crack-growth criterion, and a single value of J
cannot define a critical value. A critical toughness (e.g. J_IC or a critical crack-tip opening) would
be required in addition.

## 6. Missing or uncertain

- Assumes `J` and `C_zz` describe the same monotonic, unloaded-free state for the same material T.
- The leading-field-to-finite-body matching is only as good as the stated "negligible for r ≤ 0.2 mm".
- `T` is controlled by the ratio `J/C_zz²` (T scales roughly linearly with `(J/C_zz²)`); the angular
  coefficients depend on T only through ν_eff. Small errors in the measured C_zz propagate ~quadratically
  into T.
