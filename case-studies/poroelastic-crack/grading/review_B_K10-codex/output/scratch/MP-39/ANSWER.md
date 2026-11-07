# Suddenly Pressurized Poroelastic Crack

## Results and Limitations

The volume, pressure, final stress intensities, and loading-2 initial intensity are derived below. The true loading-1 initial tip intensity is **not completely solved**: I obtain an approximate numerical value 0.807 in the requested units, with an estimated uncertainty of about 0.01, not a rigorous error bound. Its general expression is specified by a mixed tip-layer problem rather than an evaluated closed formula. In particular, the outer undrained intensity is 1, but it must not automatically be identified with the true tip intensity. This distinction is the principal missing part of this answer.

I use Fourier-Laplace half-space solutions because they incorporate the drained crack faces and the no-flux symmetry plane outside the crack, and distinguish the inner drained tip from the outer undrained field.

| Report quantity | Value |
|---|---:|
| K1_0 | approximately 0.807, unresolved exact matching |
| K1_inf | 1 |
| V1_0 | 1.8849555921538759 |
| V1_inf | 2.5132741228718345 |
| V1_short_coeff | 1.354055000514615 |
| p_faces_0 | 0.7490636704119851 |
| K2_0 | 0 |
| K2_inf | -0.33375 |
| V2_inf | -0.8388052385084748 |

## Questions

- The prescribed face pressure is exactly zero in loading 1, so its value at the boundary is zero for all positive times. I interpret `p_faces_0` as the initial undrained pressure on the solid side, outside the vanishingly thin drainage layer, as indicated by "before it diffuses" in the question.
- The step loads have unit magnitude; reported dimensionless quantities use the normalization specified in TASK.md. Crack volume means the full displacement jump integrated from -L to L.
- No contact condition or initial aperture is specified. I use the signed linear crack solution; negative loading-2 volume denotes closure. A unilateral contact model would be a different problem.
- A finite slit maintained at nonzero pressure in infinite two-dimensional space has no bounded strictly steady harmonic pressure field vanishing at infinity. Final quantities mean limits of the transient local crack response, not a globally decaying steady pressure solution. Spatial infinity and long time do not commute.

## Reproduction

Run `python3 solve.py`. This prints all of `results.json` and regenerates that file. The analytic quantities are calculated directly; the provisional initial intensity is fitted to recorded collocation data. `python3 tip_check.py` reproduces the finest numerical collocation runs, but the coarse runs and their outputs are recorded below. Thus the initial-intensity fit is reproducible arithmetic on recorded numerical data, not a complete fresh continuum-limit computation.

## General Expressions

Write d=nu_u-nu, b=d/[alpha*(1-2*nu)], eta=alpha*(1-2*nu)/[2*(1-nu)], and tau=c*t/L^2. With dimensional K and V:

```text
K1_inf = sqrt(pi*L)
V1_0 = pi*(1-nu_u)*L^2/G
V1_inf = pi*(1-nu)*L^2/G
G*V1(t)/L^2 = pi*(1-nu_u)
             + 16*d*(1-nu_u)/[(1-nu)*sqrt(pi)]*sqrt(tau)
             + o(sqrt(tau))
p_faces_0 (outer undrained) = 2*b
K2_0 = 0
K2_inf = -eta*sqrt(pi*L)
V2_inf = -pi*alpha*(1-2*nu)*L^2/(2*G)
```

## Derivation

Initially zeta=0 outside the thin diffusion layer. The fluid-content equation gives p=-b*sigma_kk, and substitution in the strain relation replaces nu by nu_u. The elastic slit under unit opening traction has full opening

```text
w(x) = 2*(1-nu_*)/G * sqrt(L^2-x^2).
```

Its integral is pi*(1-nu_*)*L^2/G. On the faces the initial outer elastic solution has sigma_xx=sigma_yy=-1, so p=2*b. The final loading-1 pressure is zero and the compliance is drained. Although sigma_yy is negative, traction on the upper solid face is upward because its outward normal is downward.

For the transform solution let v=u_y(x,0+), zero outside the crack, so w=2v. Let T=sigma_yy, p0=p(x,0), k=abs(xi), m=sqrt(k^2+s/c), and Du=G/(1-nu_u). Decaying harmonic and diffusive modes with zero surface shear give

```text
T_hat = -Du*k*v_hat
        - 2*eta*k/(m+k)*(p0_hat-2*b*Du*k*v_hat)
p_y_hat = -m*p0_hat + 2*b*Du*k*(m-k)*v_hat.
```

At small s the first relation becomes T=-G/(1-nu)*|partial_x|v-eta*p0. Thus loading 2 is the drained elastic slit with effective opening load -eta, giving the final K and V above. At t=0+ a finite imposed pressure with zero mechanical traction has no initial singular mechanical amplitude, hence K2_0=0.

For loading 1 at large s, put ell=sqrt(c/s) and let chi be the crack indicator. The initial displacement satisfies Du*|partial_x|v0=1 on the crack. Enforcing p0=0 on the crack and p_y=0 outside gives the leading mechanical displacement correction

```text
delta(v_hat) = 4*eta*b*ell/(Du*s) * chi
delta(V_hat) = 16*eta*b*L*sqrt(c)/(Du*s^(3/2)).
```

The displacement expansion is not uniform at the endpoints, but its leading integral is the stated correction. Inversion using L^-1(s^-3/2)=2*sqrt(t/pi), and eta*b=d/[2*(1-nu)], yields the reported coefficient. It is independent of alpha when using the supplied diffusivity to normalize time.

## Unresolved Tip Matching

At every fixed positive time, arbitrarily close to the tip the mechanical singularity is drained. The outer initial field is undrained. Their intensities need not agree. In the Laplace tip layer set X=(x-L)/sqrt(c/s), r=(1-nu_u)/(1-nu), beta=2*(1-r), m=sqrt(1+k^2). Rescale the face flux to h=q/(2*b*Du), using q=p_y as the transform variable (a derivative variable, not physical Darcy flux). The exact surface multiplier is

```text
N_beta(k) = [[-k*(1-beta*k^2/[m*(m+k)]), -beta*k/[m*(m+k)]],
             [k*(m-k)/m,                 1/m]].
```

Find v,h supported on X<0 such that N_beta(|partial_X|)(v,h) vanishes on X<0, with v~sqrt(-X) as X tends to -infinity. If v~a_beta*sqrt(-X) at X=0-, then

```text
K1_0/sqrt(pi*L) = r*a_beta.
```

This is an implicit general expression, not an evaluated solution. Elastic energy matching alone cannot replace it: a contour between drained and undrained regions contains the configurational source integral of p*partial_x(zeta), which need not vanish as the layer shrinks.

## Checks and Captured Outputs

### Algebra and Normalization

An earlier independent arithmetic check compared two forms of the short-time coefficient and two forms of the steady volume, before the opening integral normalization was corrected. Its captured output was:

```text
p_faces_0 = 0.749063670411985
V1_0 = 3.76991118430775
V1_inf = 5.02654824574367
V1_short_coeff = 1.35405500051462
K2_inf = -0.33375
V2_inf = -1.67761047701695
short-coefficient reciprocity identity residual = -2.22044604925031e-16
steady-volume identity residual = 0
```

The three volumes here are **rejected**. This check illustrates its failure mode: algebraically equivalent expressions can share the same physical normalization error. The pressure, coefficient, and final K survived the corrected calculation.

A separate arithmetic comparison examined the configurational source that three candidate initial tip factors would require, without evaluating the actual source. Captured output:

```text
modulus ratio r = 0.75
sqrt(r) = 0.866025403784439
candidate K_tip/K_outer = 1: normalized hydraulic J = -0.333333333333333
candidate K_tip/K_outer = sqrt(r): normalized hydraulic J = 0
candidate K_tip/K_outer = r: normalized hydraulic J = 0.25
```

This could fail through arithmetic or an incorrect energy normalization. It does not validate any candidate, since the transition-layer source integral was not calculated.

The independent derivation compared the matrix determinant against -Du*k/m*(1-4*eta*b*k/(m+k)), compared volume formulas after integrating the full opening, and checked whether matrices at different frequencies commute. Captured output:

```text
V1_0 = 1.88495559215388
V1_inf = 2.51327412287183
V2_inf = -0.838805238508475
V1_short_coeff = 1.35405500051462
Laplace delta-volume/ell coefficient = 1.2
determinant residuals:
k=0.1 residual=-2.77555756156289e-17
k=1 residual=1.11022302462516e-16
k=10 residual=0
commutator max norm, k=0.1 and k=1 = 0.185089986588834
```

These could fail through sign, algebra, or normalization mistakes. They corrected an initial factor-of-two volume error. Determinant agreement is not verification of tip matching; the nonzero commutator only rules out constant simultaneous diagonalization.

### Collocation and Tip Resolution

I discretized the above surface relations on a periodic box of length 16L, with collocation points restricted to the crack. Unknowns are opening and p_y on the crack; opening and p_y vanish off it. I compared increasing resolution and diffusion lengths. The opening-derived tip intensity was calibrated by solving the undrained elastic problem on the same grid. Captured outputs:

```text
{"n": 512, "ell": 0.1, "K_from_opening": {"2": 0.8625018167511069, "4": 0.8602147900900093, "8": 0.8559075132085504, "16": 0.8482526835595234}, "volume": 2.0442667224357174, "relative_residual": 2.4648372002551453e-13}
{"n": 1024, "ell": 0.1, "K_from_opening": {"2": 0.8636495785864903, "4": 0.8624725077442256, "8": 0.8601884059899463, "16": 0.8558869572106955}, "volume": 2.0429251423911685, "relative_residual": 6.438592633927541e-13}
{"n": 1024, "ell": 0.05, "K_from_opening": {"2": 0.8456402932404935, "4": 0.843492940028805, "8": 0.8394496282518793, "16": 0.8322656211079595}, "volume": 1.9800714942941853, "relative_residual": 6.457811037421838e-13}
{"n": 2048, "ell": 0.05, "K_from_opening": {"2": 0.8467390386228024, "4": 0.8456332491944805, "8": 0.8434877613420592, "16": 0.8394482476543731}, "volume": 1.9794399541751084, "relative_residual": 1.0785988963958176e-12}
{"n": 2048, "ell": 0.025, "K_from_opening": {"2": 0.8337745996534414, "4": 0.8317613933948392, "8": 0.827973008281284, "16": 0.8212498650571183}, "volume": 1.9422132554242526, "relative_residual": 1.0373455207728444e-12}
{"n": 2048, "ell": 0.0125, "K_from_opening": {"2": 0.8234771072481134, "4": 0.8198934329025709, "8": 0.813541160515587, "16": 0.8034584645262108}, "volume": 1.9213397890115673, "relative_residual": 1.0773632408297936e-12}
{"n": 4096, "ell": 0.0125, "K_from_opening": {"2": 0.8254540793874435, "4": 0.8235496530717026, "8": 0.8199680396749459, "16": 0.8136191251562408}, "volume": 1.921049745437891, "relative_residual": 2.0663445842277806e-12}
```

The values tend toward the tip limit as mesh spacing decreases. A first-order mesh extrapolation, followed by a linear fit in sqrt(ell), gives 0.80698. The ell=.025 mesh correction was estimated from its nearest samples, rather than a second grid, making this particularly provisional. Periodic-image error, an incorrect extrapolation law, or inadequate tip resolution could invalidate this estimate. Tiny linear-system residuals do not establish continuum accuracy. The volume divided differences are consistent with the predicted Laplace short-time slope 1.2 but were not used as a precision validation.

The first, uncalibrated trial used the continuum opening directly at the last few grid points. Captured output:

```text
{"n": 256, "ell": 0.1, "K_from_opening": {"2": 0.9267113294444296, "4": 0.8863076378629091, "8": 0.8651506032518098, "16": 0.8470213801027662}, "volume": 2.0469454714572217, "relative_residual": 8.768221264645188e-14}
{"n": 512, "ell": 0.1, "K_from_opening": {"2": 0.9289422159062983, "4": 0.8905539543607619, "8": 0.8727649039534047, "16": 0.8592152446716217}, "volume": 2.0442667224357174, "relative_residual": 2.4648372002551453e-13}
{"n": 1024, "ell": 0.1, "K_from_opening": {"2": 0.9300890213657949, "4": 0.892804674815475, "8": 0.8770438309570403, "16": 0.8668610456943541}, "volume": 2.0429251423911685, "relative_residual": 6.438592633927541e-13}
```

This comparison could fail because endpoint collocation introduces a persistent last-grid-point bias. It motivated the same-grid elastic calibration above.

### Reproduction Output

The direct calculation prints the following numeric entries (the full printed object also contains all formula strings, the limitation, and the command):

```json
{
  "K1_0": 0.8069769386688506,
  "K1_inf": 1.0,
  "V1_0": 1.8849555921538759,
  "V1_inf": 2.5132741228718345,
  "V1_short_coeff": 1.354055000514615,
  "p_faces_0": 0.7490636704119851,
  "K2_0": 0.0,
  "K2_inf": -0.33375,
  "V2_inf": -0.8388052385084748
}
```

This checks arithmetic and JSON serialization, not the validity of the model or the unresolved initial tip limit. It could fail through a transcription error in the formulas or the recorded extrapolation inputs.

### Deliverable Consistency

I compared parsed stdout of `python3 solve.py` with the preexisting `results.json`, and checked every required key, numeric type, and formula entry. Exact captured output:

```text
required keys present: True
printed JSON equals results.json: True
numeric report entries: True
general formula coverage: True
```

This check could fail through missing keys, nonnumeric results, missing formula strings, or divergence between printed and saved values. It does not check the physics.
