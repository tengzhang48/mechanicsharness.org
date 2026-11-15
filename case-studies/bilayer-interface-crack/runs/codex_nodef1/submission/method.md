# Calculation

I used the stated plane-strain, linear-elastic model, with traction-free
crack faces and no contact. Lengths are normalized by h, and stresses by
E2, so h = 1, H = 2, E2 = 1, and E1 = 9. The two numerical solutions use
(P, M) = (1, 0) and (0, 1). No web access or external numerical packages
beyond NumPy and SciPy were used.

## Energy and amplitude

Let Ebar_i = E_i / (1 - nu_i^2). The strain energy per unit length of an
isolated layer with axial force N and centroidal moment Q is

\[
U = \frac{N^2}{2\bar E t} + \frac{6 Q^2}{\bar E t^3}.
\]

The material far ahead has zero stress. Translating the steady crack
solution releases the sum of the two remote beam energies per unit crack
advance. Thus

\[
G = \frac{P^2}{2\bar E_1 h}
  + \frac{6M^2}{\bar E_1 h^3}
  + \frac{P^2}{2\bar E_2 H}
  + \frac{6[M+P(h+H)/2]^2}{\bar E_2 H^3}.
\]

In the normalized units this becomes

\[
G(P,M)=\frac{287}{162}P^2 + 2PM + \frac{34}{27}M^2.
\]

For material 1 above the interface, define

\[
\mu_i=\frac{E_i}{2(1+\nu_i)},\qquad \kappa_i=3-4\nu_i,
\]
\[
\beta=
\frac{\mu_1(\kappa_2-1)-\mu_2(\kappa_1-1)}
     {\mu_1(\kappa_2+1)+\mu_2(\kappa_1+1)}=0.2,
\qquad
\epsilon=\frac{1}{2\pi}\log\frac{1-\beta}{1+\beta}
=-0.06453177620670408.
\]

With exactly the task's traction convention, the leading tip field is

\[
T(r)=\frac{K_h}{\sqrt{2\pi r}}(r/h)^{i\epsilon},
\qquad A=|K_h|.
\]

The plane-strain bimaterial energy relation is

\[
G=\frac{1-\beta^2}{2}
 \left(\frac{1}{\bar E_1}+\frac{1}{\bar E_2}\right)|K_h|^2
=\frac{64}{135}|K_h|^2
\]

in the chosen units. The reported normalized amplitudes therefore are
sqrt(1435/384) = 1.9331267849436742 and sqrt(85/32) =
1.6298006013006623. The energy release rates and these amplitudes use
the analytic formulas, rather than a stress extrapolation.

The phase change on increasing r tenfold is

\[
\Delta\psi=\epsilon\log(10)\frac{180}{\pi}
=-8.513574487346496\text{ degrees}.
\]

## Numerical phase calculation

`solve.py` implements displacement finite elements for the full
plane-strain elasticity equations. The truncated domain is
[-L,L] x [-2,1]. On x = -L I apply the specified remote stress as the
boundary traction. The displacement on x = L is zero, representing the
stress-free far-ahead material after removing its rigid displacement.
The outer surfaces and both crack faces have zero natural traction.
Nodes are shared across the bonded interface. The two crack faces have
distinct nodes, meeting only at the tip.

The elements are nine-node, isoparametric quadratic quadrilaterals, with
4 x 4 Gauss integration. The mesh has circular, geometrically graded
rings inside radius 0.4, transitioning to the rectangular boundary.
Angular divisions pass through the rectangle corners and material
boundaries. Collapsed quadrilaterals fill the tiny innermost ring at
the tip. The final mesh has 37,720 elements, 151,291 nodes, and 302,544
unconstrained displacement degrees of freedom. It uses L = 10, a
maximum angular interval of 2 degrees, radial growth at most 1.12,
and an innermost radius of 1e-9.

I evaluate sigma22 + i sigma12 at three interior points along each
radial interface edge, using displacement derivatives from both
materials. Their average supplies the final estimate. For each load,
I fit the complex quantity

\[
\sqrt{2\pi r}\,T(r)\exp[-i\epsilon\log(r/h)]
=K_h+c_1r+c_2r^2
\]

over 1e-4 <= r/h <= 1e-2. The constant term gives the phase at r = h of
the singular field. This extrapolation does not use the full stress
at the finite distance h. The unrounded phases from the final mesh are
52.24897442996617 and 4.287414384679429 degrees.

For M = P h, I add the complex coefficients formed from the exact
amplitudes and the fitted phases. The resulting phase is
30.436942058992894 degrees. All three reported phases are rounded to
0.01 degree.

## Checks and numerical limits

`postprocess.py` independently integrates the stated remote stress
energy density with Gauss quadrature, checks the numerical fields,
and generates `results.json` and `checks.json`. Every check passes.

| Run | L/h | Elements | Maximum angle | Radial growth | Minimum radius/h | Phase A | Phase B |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Initial | 12 | 7,296 | 5 deg | 1.22 | 1e-7 | 52.27058 | 4.40468 |
| Refined | 12 | 18,848 | 3 deg | 1.15 | 1e-8 | 52.25303 | 4.30638 |
| Shorter domain | 8 | 18,476 | 3 deg | 1.15 | 1e-8 | 52.24915 | 4.29545 |
| Final | 10 | 37,720 | 2 deg | 1.12 | 1e-9 | 52.24897 | 4.28741 |

The checks include:

- Direct integration of the remote stress energies agrees with the
  analytic energies to floating-point precision for cases A, B, and
  their sum.
- The final finite-element amplitudes differ from the exact amplitudes
  by -0.00251% and -0.02544%.
- The final upper-side and lower-side phase estimates differ by only
  0.00112 and 0.00151 degrees.
- Changing L from 8 to 12 on the refined meshes changes the two phases
  by 0.00388 and 0.01092 degrees. Refining to the final mesh changes
  them by -0.00406 and -0.01896 degrees relative to the L = 12 run.
- Moving the final fit window to [1e-5, 1e-3] changes the phases by
  -0.00242 and -0.00010 degrees.
- Fitting the logarithmic phase slope directly, without imposing its
  value, gives epsilon = -0.06451097 and -0.06453256. This independently
  checks the oscillation sign and size.
- The energy cross term requires
  cos(phase_A - phase_B) = 1/sqrt(G_A G_B), giving a phase difference
  of 47.97051474 degrees. The numerical difference is 47.96156005
  degrees, an error of -0.00895 degrees. The sign of the difference
  comes from the computed traction components.
- Superposing the fitted phases with the exact amplitudes recovers
  the combined-load energy within 0.00690%.
- The final unconstrained linear-system relative residual is
  1.96e-11. Remote load force and moment resultants cancel to about
  1e-15 in the normalized units.

The finite domain, finite mesh, and near-tip extrapolation remain
numerical approximations. The convergence checks support phase errors
of a few hundredths of a degree, comfortably below the requested
0.5-degree tolerance; this is not a rigorous error bound. I have not
compared with an independent reference solution. No contact model was
introduced, as specified in the task.

## Reproduction

Run these commands from this folder. Each solve writes its fit data
to JSON and its sampled interface stresses to a matching NPZ file.

```sh
OPENBLAS_NUM_THREADS=1 python3 solve.py --output fem_initial.json
OPENBLAS_NUM_THREADS=1 python3 solve.py --length 12 --angle 3 --ratio 1.15 --rmin 1e-8 --output fem_refined.json
OPENBLAS_NUM_THREADS=1 python3 solve.py --length 8 --angle 3 --ratio 1.15 --rmin 1e-8 --output fem_length8.json
OPENBLAS_NUM_THREADS=1 python3 solve.py --length 10 --angle 2 --ratio 1.12 --rmin 1e-9 --output fem_fine.json
OPENBLAS_NUM_THREADS=1 python3 postprocess.py
```

All code, saved numerical data, checks, and requested deliverables are
contained in this folder.
