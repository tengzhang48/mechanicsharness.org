# Derivation notes

## Beam-state stresses

Layer 1 (`0 < y < h`):

```
s11_1(y) = -P/h - 12 M (y - h/2)/h^3
```

Layer 2 (`-H < y < 0`):

```
s11_2(y) = P/H + 12 (M + P (h+H)/2) (y + H/2)/H^3
```

Both layers have `s12 = s22 = 0` far behind.  Direct integration gives net
axial force zero and net moment zero about the origin, and the far-ahead state
is therefore stress free.

## J integral

For a contour enclosing the tip, running along the free surfaces and the
cross-sections far behind and ahead, the only non-vanishing contributions are
the two far-behind vertical cuts:

```
G = int_0^h s11_1^2/(2 E1*) dy + int_{-H}^0 s11_2^2/(2 E2*) dy
```

with `E_i* = E_i/(1-nu^2)`.  The far-ahead cut is stress free.  This gives

```
G_A = P^2 [ 1/(18h E1*) + 1/(2H E2*) + 3(h+H)^2/(H^3 E2*) ]
G_B = 6 M^2/(E1* h^3) + 6 M^2/(E2* H^3)
```

For `E1 = 9 E2`, `h/H = 1/2`, `nu = 1/3` and `E2 = H = 1`:

```
G_A = 8.28125 P^2,   G_A E2 h/P^2 = 4.140625
G_B = 12.75 M^2,     G_B E2 h^3/M^2 = 1.59375
```

## Interface-crack amplitude

The standard relation for the complex stress intensity factor
`T(r) = K (2 pi r)^(-1/2) (r/h)^(i eps)` is

```
G = C_G |K|^2,
C_G = [(kappa+1)/mu1 + (kappa+1)/mu2]/(16 cosh^2(pi eps)).
```

With `kappa = 5/3`, `mu1 = 27/8`, `mu2 = 3/8`,
`cosh(pi eps) = 5/(2 sqrt(6))`, this is exactly `C_G = 64/135`.  Hence

```
case A: |K| = sqrt(8.28125 * 135/64) = 4.17950496097
        A sqrt(h)/P = 4.17950496097 * sqrt(1/2) = 2.9553562999
case B: |K| = sqrt(12.75 * 135/64) = 5.18599375723
        A h^(3/2)/M = 5.18599375723 * (1/2)^(3/2) = 1.8335256765
```

## Oscillatory exponent

Dundurs parameters for equal Poisson ratios:

```
alpha = (E1* - E2*)/(E1* + E2*) = 0.8
beta  = (1/2) (1-2 nu)/(1-nu) (mu1-mu2)/(mu1+mu2) = 0.2
eps   = (1/(2 pi)) ln[(1-beta)/(1+beta)] = -0.0645317762067
```

The Airy eigenfunction was checked independently: for
`f(theta) = A exp(i s theta) + B exp(-i s theta)
            + C exp(i(s+2)theta) + D exp(-i(s+2)theta)`,
the 8x8 system enforcing traction-free faces at `theta = +/- pi` and
continuity of traction and displacement at `theta = 0` has determinant zero
at `s = -1/2 +/- i eps`.

The phase change per decade is `eps ln 10 = -8.513574487 deg`.

## Numerical calibration

A manufactured problem was solved with the same FE code: the exact analytic
interface-crack eigenfunction was imposed as tractions on the boundary of a
finite box.  Comparing the FE crack-opening displacement with the exact
`K` of that eigenfunction calibrates the complex constant in

```
delta_y + i delta_x = C K r^(1/2+i eps).
```

The calibrated phase is `arg C = -arctan(2 eps) = +7.36 deg`, i.e. the
standard factor `1/(1+2 i eps)` in the crack-opening formula.  This
calibration was used for the jump-based phase extraction.
