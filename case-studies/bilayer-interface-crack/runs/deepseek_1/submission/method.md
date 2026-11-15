# Method and results

## 1. Units, geometry and material data

Lengths are measured in units of the upper-layer thickness, so

```
h = 1,   H = 2,   E2 = 1,   E1 = 9,   nu1 = nu2 = 1/3.
```

Plane strain is used throughout.  The shear moduli and Kolosov constants are

```
mu_k  = E_k / (2(1+nu_k))          =>  mu1 = 27/8,  mu2 = 3/8
kappa = 3 - 4 nu                   =>  kappa1 = kappa2 = 5/3.
```

The Dundurs parameters are

```
beta  = (mu1(kappa2-1) - mu2(kappa1-1)) / (mu1(kappa2+1) + mu2(kappa1+1)) = 0.2
alpha = (E1* - E2*)/(E1* + E2*) = 0.8,     E_k* = E_k/(1-nu_k^2).
```

The bimaterial constant requested in the task is

```
eps = (1/2pi) ln[ (kappa1/mu1 + 1/mu2) / (kappa2/mu2 + 1/mu1) ]
    = (1/2pi) ln(2/3)
    = -0.0645317762067041.
```

This is exactly the value used in the numerical work below (it also satisfies
`tanh(pi*eps) = -beta`).

## 2. Energy release rate from beam theory / the J-integral

Far behind the tip (`x1 -> -inf`) each layer is a separate beam with the
prescribed stress.  That stress state is an exact elasticity solution for the
two separated layers.  Far ahead (`x1 -> +inf`) the bonded bilayer is
stress-free.  Take a J-integral contour which follows the two free surfaces,
the two traction-free crack faces and vertical cuts far behind and far ahead of
the tip.  The only non-zero contribution comes from the far-behind cut, where
`u_{i,x}` is the beam strain.  The result is simply the strain energy per unit
length of the two beams far behind:

```
G = int_{-H}^{0} sigma11^2/(2 E2*) dx2 + int_{0}^{h} sigma11^2/(2 E1*) dx2.
```

For plane strain the uniaxial compliance is `1/E_k*`.  With the stress
distributions given in the problem statement, the cross terms between axial
force and bending moment integrate to zero, and

```
G_A = P^2/(2 E1* h) + P^2/(2 E2* H)
      + 6 [P(h+H)/2]^2/(E2* H^3)
    = 287/162 * P^2/(E2 h)                 (M = 0)

G_B = 6 M^2/(E1* h^3) + 6 M^2/(E2* H^3)
    = 34/27 * M^2/(E2 h^3)                 (P = 0).
```

Numerically,

```
G_A E2 h/P^2 = 1.7716049382716048
G_B E2 h^3/M^2 = 1.2592592592592593.
```

This is exact for the semi-infinite crack in the strip; it does not rely on any
numerical approximation.  As an independent check, the same values were
recovered from the finite-element solution by evaluating the J-integral on
vertical cuts behind the tip (see section 4.3): the agreement was at the level of
`1e-8` relative.

## 3. Magnitude of the complex stress intensity factor

For an interface crack the standard relation between the energy release rate
and the complex stress intensity factor defined by

```
sigma22 + i sigma12 = K r^(i eps) / sqrt(2 pi r)     (r -> 0)
```

is

```
G = (1 - beta^2) |K|^2 / E*,       1/E* = (1/2)(1/E1* + 1/E2*).
```

Here `E* = 2.025 E2` and `1 - beta^2 = 0.96`.  Therefore

```
|K_A| sqrt(h)/P   = sqrt( (287/162) * (135/64) ) = 1.933126784943674
|K_B| h^(3/2)/M   = sqrt( (34/27)  * (135/64) ) = 1.629800601300662.
```

The factor `135/64 = 2.025/0.96` is `E*/(1-beta^2)`.  The relation was checked
directly against the exact near-tip eigenfield: for a unit `K` the J-integral of
the eigenfield is `0.4740747`, in agreement with `(1-beta^2)/E* = 0.4740741` to
about `1e-6` relative.

## 4. Phase angles

### 4.1 Finite-element model

A plane-strain displacement finite-element program was written for this problem
(`fem_bilayer.py`).  It uses high-order `Q_p` isoparametric quadrilaterals
(p = 2 and 3) on a structured, geometrically graded mesh:

* domain `-L <= x1 <= L`, `-H <= x2 <= h`, with `L = 10 h` (checks with `L = 12 h`);
* layer 1 occupies `0 < x2 < h`, layer 2 `-H < x2 < 0`;
* the crack is `x1 < 0, x2 = 0`: the upper and lower crack-face nodes are
  duplicated and traction-free;
* the interface `x1 > 0, x2 = 0` is perfectly bonded (shared nodes);
* the outer surfaces `x2 = h` and `x2 = -H` are traction-free;
* at `x1 = -L` the exact beam tractions `t = -sigma11(y) e_x` are applied
  (consistent nodal forces); at `x1 = +L` the boundary is traction-free;
* three displacement constraints at the far right remove the rigid-body modes;
* the mesh is geometrically refined towards the tip in both directions; the
  finest runs had about `2.9e5` degrees of freedom, and the sparse system was
  solved directly.

The element stiffness was checked against a uniform-strain patch, and a
uniform-tension FEM test reproduced the exact uniform stress to better than
`1e-5`.

### 4.2 Crack-opening-displacement (COD) extraction

Let `r = -x1 > 0` be the distance behind the tip and let

```
delta2 + i delta1 = (u2^upper - u2^lower) + i (u1^upper - u1^lower).
```

For the interface crack the asymptotic COD is

```
delta2 + i delta1 = [ 8 / (E* cosh(pi eps)) ] [ 1/(1 + 2 i eps) ]
                    sqrt(r/(2 pi)) K r^(i eps) + O(r^(3/2+i eps)).
```

The `1/(1+2 i eps)` factor is important for the phase.  It was determined by
constructing the exact Airy eigenfunction of the interface-crack problem
(`asymptotic.py`; eigenvalue `lambda = 1/2 + i eps`), normalising it so that its
stress combination at `theta = 0` is exactly `K r^(i eps)/sqrt(2 pi r)`, and
evaluating its displacement jump.  The ratio of the exact jump to the
uncorrected formula is exactly `1/(1+2 i eps)` (checked to `1e-6` over several
decades in `r` and several complex `K`).  For `eps = 0` the factor is one and the
standard homogeneous result is recovered.  In the homogeneous limit the
eigenfield reproduces the standard mode-I angular stress field (to the accuracy
of the degenerate homogeneous eigenvector); for the bimaterial case its
J-integral equals `(1-beta^2)|K|^2/E*` to about `1e-6`, which validates both its
stress and displacement parts.

From the FEM solution, the nodal COD values were used to form

```
Q(r) = [ delta2 + i delta1 ] sqrt(2 pi/r) r^(-i eps) / C_eff,
C_eff = 8/[E* cosh(pi eps)(1 + 2 i eps)].
```

`Q(r)` tends to `K` as `r -> 0`; the leading correction is `O(r)`.  A
low-degree polynomial in `r` was fitted over the window `r/h in [0.01, 0.1]`
(and `[0.02, 0.1]`) and extrapolated to `r = 0`.  The windows and polynomial
degree were varied; the spread of the extracted phase is below `0.05 deg`.
The same extraction was repeated on several meshes.  For cases A and B the
fitted magnitudes are `1.928` and `1.627`, i.e. about `0.2-0.3%` below the
closed-form values `1.9331` and `1.6298`; this is the expected small
discretisation/extrapolation bias, so the reported magnitudes are taken from the
closed-form `G`-relation rather than from the fit.

### 4.3 Results of the phase extraction

Central COD fits (finest mesh, average of the two central windows):

```
case A (P>0, M=0):   psi = 52.27 deg
case B (P=0, M>0):   psi =  4.25 deg
```

Mesh convergence of the extracted phases (COD fit, `r` window `[0.01,0.1] h`):

| mesh | dofs | psi_A (deg) | psi_B (deg) |
|---|---|---|---|
| p=2, nx=100, ny=30/36 | 1.07e5 | 52.31 | 4.16 |
| p=3, nx=60,  ny=24/30 | 1.18e5 | 52.30 | 4.20 |
| p=3, nx=80,  ny=30/36 | 1.92e5 | 52.28 | 4.24 |
| p=3, nx=100, ny=36/44 | 2.90e5 | 52.27 | 4.25 |

Independent cross-checks on the finest meshes:

1. **Stress extrapolation.**  The stresses `sigma22 + i sigma12` on the bonded
   interface ahead of the tip were used directly in the definition
   `K = (sigma22 + i sigma12) sqrt(2 pi r) r^(-i eps)` and extrapolated to
   `r -> 0`.  This gives `psi_A = 52.20 - 52.28 deg`, `psi_B = 4.13 - 4.22 deg`.
2. **Interaction (M) integral.**  A domain interaction integral with the exact
   Airy eigenfield as the auxiliary field gave, on a coarse mesh,
   `psi_A = 52.06 - 52.40 deg`, `psi_B = 4.12 - 4.40 deg` and `|K|` within about
   `0.3%` of the exact values.  The variation with the integration domain is the
   numerical uncertainty of this check.
3. **Homogeneous limit.**  Setting `E1 = E2`, `h = H` and using load case B
   gives a symmetric mode-I problem whose exact solution has
   `K = sqrt(12) = 3.46410` and `psi = 0`.  The same FEM/COD procedure returned
   `psi = 0` (to `1e-8 deg`) and `|K| = 3.456 - 3.459` (0.2% low, the expected
   discretisation bias of the extrapolation).  This validates the COD constant
   and the phase convention.
4. **Energy check.**  The FEM J-integral on vertical cuts behind the tip gave
   `J_A = 1.77160493`, `J_B = 1.25925926` and, for the homogeneous symmetric
   test, `J = 10.66666668`, in agreement with the exact beam energies to about
   `1e-8` relative.

### 4.4 Both loads with `M = P h`

Because the problem is linear, `K_total = K_A(P) + K_B(M)` with `M = P h`.
Using the exact magnitudes and the fitted phases gives

```
psi_both = arg( K_A + K_B ) = 30.43 deg.
```

A direct FEM run with `P = 1, M = 1` (i.e. `M = P h`) gave, from the same COD
fit, `psi = 30.42 deg` (stress extrapolation `30.35 - 30.42 deg`), in agreement
with the superposition.  The reported value is `30.42 deg`.

## 5. Reported numbers

`results.json` contains

```
eps                 = -0.0645317762067041
case A:
  G E2 h / P^2      = 1.7716049382716048
  |K| sqrt(h) / P   = 1.933126784943674
  psi               = 52.27 deg
case B:
  G E2 h^3 / M^2    = 1.2592592592592593
  |K| h^(3/2) / M   = 1.629800601300662
  psi               = 4.25 deg
psi_both            = 30.42 deg
```

## 6. Assumptions, checks and remaining uncertainty

* Plane strain, small strains, no body forces, no crack-face contact, perfect
  bonding for `x1 > 0`, traction-free outer surfaces and crack faces.
* The far-field beam tractions are applied at a finite distance `L = 10 h`
  (`12 h` in one check).  The J-integral is independent of the vertical cut
  position to `1e-8`, confirming that the boundary-layer/truncation effects are
  negligible for the energy and the phase.
* `eps`, the two energy release rates and the two `|K|` values are obtained in
  closed form and are not subject to numerical error beyond the validation
  checks described above.
* The phase angles are numerical.  The spread between fit windows, meshes and
  independent extraction methods is at most about `0.2 deg`; the quoted values
  are estimated to be accurate to about `0.1 deg`, well inside the required
  `0.5 deg`.  There is no closed-form reference for the present geometry and
  Dundurs parameters used here, so this small numerical uncertainty is the main
  remaining unverified item.
* No crack-face contact or friction was included, as specified.

## 7. Code

* `fem_bilayer.py` — mesh generation, high-order plane-strain FEM, beam
  tractions, constraints, COD/stress/J-integral/interaction-integral
  post-processing.
* `asymptotic.py` — exact Airy eigenfunction of the interface crack, complex
  `K` normalisation, `1/(1+2 i eps)` COD factor, auxiliary fields used in the
  interaction-integral checks.
* `analyze_case.py` — solves a load case and performs the COD/stress fits.
* `converge.py` — mesh-convergence runs.
* `finalize_results.py` — evaluates `eps`, `G`, `|K|` in closed form, fits the
  saved COD data and writes `results.json`.
* `data/cod_A.npz`, `data/cod_B.npz`, `data/cod_C.npz` — saved COD data for
  cases A, B and the combined load `P=1, M=1` on the finest mesh.

To regenerate the saved data and the final JSON:

```
python3 analyze_case.py A     # ~40 s, writes data/cod_A.npz
python3 analyze_case.py B     # ~40 s, writes data/cod_B.npz
python3 analyze_case.py C     # ~40 s, writes data/cod_C.npz (P=1, M=1)
python3 finalize_results.py   # writes results.json
```
