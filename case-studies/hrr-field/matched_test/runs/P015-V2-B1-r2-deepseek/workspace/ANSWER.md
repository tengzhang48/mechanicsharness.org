# Inference of the stiffness parameter T from a transverse strain field

Units: mm, N, MPa. Strains are tensor components (engineering shear = 2 x tensor shear).

## 1. Requested numbers

| Quantity | Value |
|---|---|
| `T` | **4816.47 MPa** (4.8165e3 MPa) |
| `C_rr(60°) = lim sqrt(r) eps_rr` | **3.9420e-3 mm^1/2** |
| `C_θθ(60°) = lim sqrt(r) eps_θθ` | **6.3318e-4 mm^1/2** |
| `C_γrθ(60°) = lim sqrt(r) 2 eps_rθ` | **5.7311e-3 mm^1/2** |
| `C_y(90°) = lim u_y / sqrt(r)` | **9.1390e-3 mm^1/2** |
| Screening interval on θ = 0 | **[0.010000, 0.061325] mm** (non-empty) |

Auxiliary quantities used: effective modulus `E_eff = T = 4816.47 MPa`, effective Poisson ratio
`ν_eff = 0.486239`, effective shear modulus `G_eff = 1620.36 MPa`, stress-angular amplitude
`K_A = 24.7639 MPa·mm^1/2` (defined by `σ_ij = r^(-1/2) Σ_ij(θ)`, so `Σ_yy(0)=K_A`).

## 2. Model chosen and why

**Leading field.** The crack is a straight, stationary, traction-free slit (β = π), mode I by the
symmetric tension. Close to the tip the effective constitutive response is *linear* (see below),
so the leading in-plane field is the Williams mode-I field, and the exponent is fixed by J:
`p = 1/2`, i.e. `σ ~ r^(-1/2)`, `ε ~ r^(-1/2)`, `u ~ r^(1/2)` — exactly the `√r` coefficients
requested. With `σ_zz = σ_xz = σ_yz = 0` this is **plane stress**, and the in-plane stress angular
functions are universal (independent of the elastic constants):

```
Σ_rr = K_A [ 5/4 cos(θ/2) - 1/4 cos(3θ/2) ]
Σ_θθ = K_A [ 3/4 cos(θ/2) + 1/4 cos(3θ/2) ]
Σ_rθ = K_A [ 1/4 sin(θ/2) + 1/4 sin(3θ/2) ]        σ_ij = r^(-1/2) Σ_ij(θ)
```

**Constitutive reduction (as r → 0).** Write `S = σ_kk`, `s_ij` the deviator. The given total-strain
law is

```
ε_ij = A_hyd S δ_ij + A_dev s_ij + (3/2)(1/T - 1/E) ( <1 - σ_y/σ_e> - 1 ) s_ij
A_hyd = (1-2ν)/(3E),      A_dev = 3/(2T) + (ν-1/2)/E .
```

Because the stress is singular, `σ_e = O(r^(-1/2)) → ∞`, and the last term is of lower order
(`<1-σ_y/σ_e>-1 = -σ_y/σ_e = O(r^(1/2))` relative to the `s` term). Hence to leading order

```
ε_ij = A_hyd S δ_ij + A_dev s_ij .
```

This is an **isotropic linear-elastic law** with `1/(2G_eff) = A_dev` and `1/(9K) = A_hyd`:

```
E_eff = T                       (exact: uniaxial slope 1/E_eff = A_hyd + (2/3)A_dev = 1/T)
ν_eff = 1/2 + T(ν - 1/2)/E = 0.5 - 0.2 T/E .
```

The result `E_eff = T` is not an approximation: for any stress with `σ_e > σ_y`, the uniaxial
response of the full law is `ε = σ/E + (1/T-1/E)(σ-σ_y)`, i.e. exactly linear with tangent modulus
`T`, so `T` *is* the asymptotic Young modulus. The yield term therefore only corrects *subleading*
(regular) terms; it cannot change the `r^(-1/2)` coefficients.

**Treatment of ε_zz.** σ_zz = 0 (plane stress), but ε_zz is free and non-zero. It follows from the
reduced law: `ε_zz = A_hyd S - (A_dev/3) S = (A_hyd - A_dev/3) S`. In effective-elastic language,
`ε_zz = -(ν_eff/E_eff)(σ_xx+σ_yy) = -(ν_eff/E_eff) S`, which is what makes `C_zz(0) ≠ 0`.

**Amplitude from J.** For a linear-elastic material in plane stress, `J = K_std^2/E_eff` with the
standard SIF `K_std = sqrt(2π) K_A`, so

```
J = 2π K_A^2 / E_eff      =>      K_A = sqrt( J E_eff / (2π) ) = 24.7639 MPa·mm^1/2  (for T above).
```

## 3. Determination of T

At θ = 0, `σ_xx = σ_yy = K_A r^(-1/2)`, `σ_e = K_A r^(-1/2)`, and

```
C_zz(0) = -2 (ν_eff/E_eff) K_A = -2 ν_eff sqrt( J / (2π E_eff) ) .
```

Setting `C_zz(0) = -0.005000`, `J = 0.800`, `E_eff = T`, `ν_eff = 0.5 - 0.2 T/E`:

```
(0.5 - 0.2 T/70000) / sqrt(T) = 0.0025 sqrt(2π/0.800)  = 7.00624e-3
=>  0.5 - T/350000 = 7.00624e-3 sqrt(T)  =>  T = 4816.472 MPa .
```

## 4. θ = 60° coefficients

At θ = 60° (`cos(θ/2)=√3/2`, `cos(3θ/2)=0`, `sin(θ/2)=1/2`, `sin(3θ/2)=1`):

```
Σ_rr = 1.082532 K_A,  Σ_θθ = 0.649519 K_A,  Σ_rθ = 0.375 K_A,  S = 1.732051 K_A .
```

With `ε_rr = (Σ_rr - ν_eff Σ_θθ)/E_eff`, `ε_θθ = (Σ_θθ - ν_eff Σ_rr)/E_eff`,
`ε_rθ = (1+ν_eff) Σ_rθ/E_eff` (plane stress, tensor shear):

```
C_rr   = (Σ_rr - ν_eff Σ_θθ)/E_eff = 3.9420e-3 mm^1/2
C_θθ   = (Σ_θθ - ν_eff Σ_rr)/E_eff = 6.3318e-4 mm^1/2
C_γrθ  = 2(1+ν_eff) Σ_rθ/E_eff     = 5.7311e-3 mm^1/2
```

## 5. C_y(90°)

At θ = 90°, `e_r = (0,1)`, so `u_y = u_r`. The mode-I angular displacement (plane stress,
`κ_eff = (3-ν_eff)/(1+ν_eff) = 1.69137`) is

```
u_r(r,θ) = (K_A/(2G_eff)) sqrt(r) (κ_eff - cos θ) cos(θ/2) ,
```

so

```
C_y(90°) = (K_A/(2G_eff)) κ_eff / sqrt(2) = 9.1390e-3 mm^1/2 .
```

(Equivalent route: `du_r/dr = ε_rr` on the ray gives `u_r/sqrt(r) = 2 ε_rr(90°)`, which reproduces
the same number.)

## 6. Screening interval on θ = 0

At θ = 0 the three principal strains are `ε_xx = ε_yy = (1-ν_eff)σ/E_eff` and
`ε_zz = -2ν_eff σ/E_eff`, with `σ = K_A/sqrt(r)`.

* `σ_e ≥ 10 σ_y = 100 MPa`:  `K_A/sqrt(r) ≥ 100`  →  `r ≤ (K_A/100)^2 = 0.0613252 mm`.
* `max|principal strain| ≤ 0.05`: the largest is `|ε_zz| = (2ν_eff σ/E_eff) = |C_zz(0)|/sqrt(r) = 0.005/sqrt(r)`,
  so `0.005/sqrt(r) ≤ 0.05` → `r ≥ 0.0100000 mm`. (Equivalently `σ ≤ 247.639 MPa`.)
* `r ≤ 0.200 mm`.

Intersection: **`0.010000 ≤ r ≤ 0.061325` mm**, non-empty.

## 7. Verification (quantitative, not the C_zz substitution)

1. **Direct J integral.** I built the physical field (`K_A` above, effective material `E_eff=T`,
   `ν_eff`), computed the Cartesian displacement and stresses including the polar-basis
   derivatives, and evaluated the problem's `J = ∫(W N_x - t_i u_{i,x}) ds` on circular arcs:
   `J = 0.8000000000` at R = 0.05, 0.10, 0.20 and 1.00 mm (independently at four radii). This
   confirms the amplitude/K–J normalization.
2. **Topic constraint checks** (`knowledge/checks/crack_tip_checks.py`, `check_field`, effective
   linear law with `ν_eff`, `p = 1/2`): all PASS — faces traction-free (1e-23), mode-I symmetry
   (2e-16), equilibrium angular (1e-10) and Cartesian FD (1e-8), compatibility (St. Venant 4e-7;
   displacements from strains 2e-12), finite energy, crack exponent fixed by J, crack faces open.
3. **Full nonlinear law check.** Substituting the leading stress field into the *original* law and
   letting `σ_y/σ_e → 0` recovers the reduced coefficients: at θ=0, `sqrt(r) ε_zz` = -0.0048067,
   -0.0049389, -0.0049807, -0.0049939 at r = 1e-2, 1e-3, 1e-4, 1e-5 (→ -0.005, bracket factor
   `1-σ_y/σ_e` = 0.960, 0.987, 0.996, 0.999); at θ=60° the three coefficients similarly converge
   to 3.9420e-3, 6.3318e-4, 5.7311e-3. This verifies the constitutive reduction, including the
   claim that the yield term is subleading.

## 8. What the screens establish, and crack initiation

The three screens select the annular (here: radial) window on θ=0 where the leading singular
field is simultaneously (i) deep inside the plastic zone (`σ_e ≥ 10σ_y`, so the asymptotic linear
response dominates), (ii) still small-strain (`max|principal strain| ≤ 0.05`), and (iii) inside
the region (`r ≤ 0.2 mm`) where matching to the outer body is stated to be negligible. They are
consistency/applicability checks on the asymptotic field, **not** error estimates for the full
finite-r solution.

The supplied data do **not** determine crack initiation. J and T (through C_zz) fix the strength
(field amplitude) and the material's asymptotic stiffness, but no fracture criterion (critical J,
critical crack-opening, toughness, or a process-zone/damage model — explicitly excluded) is given.
Per the topic (§6 item 9), the amplitude is a field parameter, not a fracture criterion; a notch/
crack onset condition would require an additional energy/strength criterion.

## 9. Assumptions changed and why

* **Plane strain → plane stress.** The problem states `σ_zz = σ_xz = σ_yz = 0`, so plane stress.
* **Using the asymptotic (large-σ_e) branch of the law as the leading field.** Justified because
  the stress is singular; the yield cutoff contributes only regular terms (verified in §7.3).
* No material parameter other than T is used; extension/tension is mode I by symmetry.

## 10. Missing / uncertain

* T is inferred from a single measured coefficient plus J; the two are mutually consistent by
  construction, so the result is only as good as `C_zz(0) = -0.005000` and `J = 0.800` (given as
  exact "reports"). With the stated 1% tolerance on T this is not a concern for the answer.
* The interval endpoints are exact for the leading field; the upper end lies below 0.2 mm, so the
  `r ≤ 0.2` condition is not active.
* Displacements are fixed only up to the stated absence of additive rigid translation; a rotational
  rigid mode is excluded by mode-I symmetry.
