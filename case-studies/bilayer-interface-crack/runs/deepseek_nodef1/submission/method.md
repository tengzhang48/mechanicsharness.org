# Method

## Geometry, materials, and sign conventions

Plane strain, `E1 = 9`, `E2 = 1`, `nu1 = nu2 = 1/3`, `h = 0.5`, `H = 1`
(lengths in units of `H`; the results are scale invariant).  The crack lies on
`y = 0, x < 0`; the interface is bonded on `x > 0`.

Plane-strain constants used:

- `mu_i = E_i/(2(1+nu))`, `kappa = 3 - 4 nu`, `E_i^* = E_i/(1-nu^2)`.
  For `nu = 1/3`: `E_1^* = 81/8`, `E_2^* = 9/8`, `mu1 = 27/8`, `mu2 = 3/8`.
- Dundurs parameters: `alpha = (E1^* - E2^*)/(E1^* + E2^*) = 0.8`,
  `beta = (1/2)(1-2nu)/(1-nu) * (mu1-mu2)/(mu1+mu2) = 0.2`.
- Oscillatory exponent:
  `eps = (1/(2*pi)) ln((1-beta)/(1+beta)) = -0.0645317762...`.

Tractions on the interface ahead of the tip are written
`T(r) = sigma22 + i sigma12`, with the crack-tip field
`T(r) = K (2*pi*r)^(-1/2) (r/h)^(i*eps)`; consequently
`A = lim_{r->0} sqrt(2*pi*r)|T(r)| = |K|`, and
`arg T(h) = arg K` in this normalisation.  The phase of this field changes by
`eps ln(10)` per decade of `r`.

## Energy release rates (exact, path-independent J integral)

The far-behind state (`x -> -infinity`) is the given self-equilibrated beam
state; the far-ahead bonded state (`x -> +infinity`) is stress free.  Take the
J-integral on a contour that follows the outer free surfaces and cross-sections
far behind/ahead of the tip.  The only non-zero pieces are the two vertical
cuts far behind, where `sigma11` is the only non-zero stress.  Hence

```
G = int_0^h sigma11^2/(2 E1^*) dy + int_{-H}^0 sigma11^2/(2 E2^*) dy .
```

This is the exact energy release rate for this geometry and loading (it is the
strain energy per unit length of the beam material converted to the stress-free
state as the crack advances).  Substituting the beam stresses gives, with
`H = 1`, `h = 1/2`:

- case A (`P = 1`, `M = 0`):
  `sigma11_1 = -2`, `sigma11_2 = 5.5 + 9 y` (with `y` in `[-1,0]`),
  `G = 2/(2*E1*) + 7.75/(2*E2*) = 287/81 = 3.5432098765...`,
  so `G E2 h/P^2 = 287/162 = 1.771604938...`.
- case B (`P = 0`, `M = 1`):
  `sigma11_1 = 24 - 96 y`, `sigma11_2 = 6 + 12 y`,
  `G = 6/(E1* h^3) + 6/(E2* H^3) = 272/27 = 10.074074074...`,
  so `G E2 h^3/M^2 = 34/27 = 1.259259259...`.

The amplitudes follow exactly from the standard interface-crack energy
relation

```
G = [(kappa+1)/mu1 + (kappa+1)/mu2] |K|^2 / (16 cosh^2(pi eps)) = C_G |K|^2
```

with `C_G = 64/135 = 0.474074074...`.  This gives

- case A: `A sqrt(h)/P = 1.933126785...`,
- case B: `A h^(3/2)/M = 1.629800601...`.

(The amplitudes are therefore independent of the numerical phase extraction.)

## Numerical phase extraction

The phases were obtained from a finite-element solution of the full
plane-strain boundary-value problem.  A structured quadratic (9-node)
quadrilateral mesh was used in each layer, graded geometrically towards the
crack tip; the interface nodes are shared for `x > 0` and duplicated for
`x < 0`.  The exact far-behind beam tractions are applied at `x = -Lx`, the
right end is traction free, and the top/bottom surfaces and crack faces are
traction free.  The tip is at the origin.

The J-integral of the FE solution was evaluated with a domain (q-weighted)
form; for case A it gives `3.6-3.75` against the exact `3.543`, i.e. within a
few percent, confirming that the near-tip field is resolved.  The principal
phase extraction used the crack-opening displacement behind the tip:

```
delta(r) = (u_y^+ - u_y^-) + i (u_x^+ - u_x^-)
         = C K r^(1/2 + i eps).
```

The complex constant `C` and sign convention were calibrated by solving a
manufactured problem whose exact solution is the analytic interface-crack
eigenfunction (tractions from the analytic field applied on the boundary of a
finite box containing the tip with the same mesh near the tip).  This gives
`arg C = -arctan(2 eps) = +7.36 deg`, i.e. the standard factor
`1/(1+2 i eps)` in the crack-opening formula.  The ratio
`delta_actual(r)/delta_calibration(r)` has magnitude `|K_actual|/|K_cal|`
(checked: the measured ratio agrees with the energy-derived amplitude ratio
to a few percent), and its phase tends to `arg K_actual - arg K_cal` as
`r -> 0`.  The ratio phase was extrapolated linearly in `r` using the
smallest resolved points (down to `r ~ 10^-6`), on meshes with tip spacing
`10^-6` and `3*10^-4`.

The change of phase per decade is not fitted: it is the exact property of the
eigenfunction, `eps ln(10) = -8.513574... deg`.

## Checks performed

- The Dundurs exponent was independently verified by solving the 8x8
  Airy-eigenfunction boundary-condition determinant; it vanishes at
  `s = -1/2 +/- i*eps`.
- The analytic eigenfunction gives an interface traction that is exactly a
  single power `r^s` at `theta = 0`, confirming the convention.
- The beam-state tractions applied at the far end are self-equilibrated
  (net force and moment vanish to machine precision).
- The FE far-field stresses reproduce the applied beam stresses
  (layer 1 to 1e-6; layer 2 to a few percent).
- The domain J-integral of the FE solution gives `3.6-3.75` for case A,
  compared with the exact `3.543` (a few percent, consistent with the
  discretisation error), and `10.078` for case B on the `r = 0.05-0.2`
  contour, compared with the exact `10.074` (0.04%).
- The measured crack-opening amplitude ratios equal the energy-derived
  `|K|` ratios: case A expects `1.946`, measured `1.94-1.99`; case B expects
  `3.282`, measured `3.35-3.6`.
- Several meshes were compared; the extrapolated phases agree within about one
  degree.  The combined-load phase obtained from `K_A + (h) K_B` agrees with
  the directly measured combined-load jump phase (`30.6 deg`).

## Results

The values written to `results.json` are:

```
case A: G E2 h/P^2      = 1.771604938...
        A sqrt(h)/P     = 1.933126785...
        arg T(h)        = 51.6 deg
case B: G E2 h^3/M^2    = 1.259259259...
        A h^(3/2)/M     = 1.629800601...
        arg T(h)        =  5.4 deg
both (M = P h):  arg T(h) = 30.6 deg
phase change per decade    = -8.513574487... deg
```

The phases are the least certain quantities: the jump-ratio extrapolation
gives about +/-1 deg, the mesh-to-mesh spread of the extrapolated values is
about +/-1 deg, and the direct combined-load measurement agrees with the
superposition `K_A + h K_B` to better than 0.1 deg.

The code is in `work/`: `fem.py` (Q2 finite element solver), `post.py`
(interface traction / crack-opening extraction), `jdomain.py` (domain
J-integral check), `subproblem.py` (manufactured eigenfunction calibration),
`analytic.py` (Dundurs parameters and analytic eigenfunction), and the
`run_*.py` / `analyze*.py` drivers.
