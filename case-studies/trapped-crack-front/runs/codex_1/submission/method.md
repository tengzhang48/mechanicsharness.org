# First-order solution

All calculations use the specified infinite-solid loading, a crack confined to
the plane, and the local normal and tangent in the question. Terms beyond first
order in front amplitude or in epsilon are omitted. No web sources or external
data were used. Run `python3 solve.py` in this folder to reproduce `results.json`
and all the numerical checks.

## 1. Elastic response of a planar perturbed front

Let `h(z) = a(z) - a_bar`, and define the line operator

\[
\mathcal J[h](z)=\frac{1}{2\pi}\operatorname{PV}
 \int_{-\infty}^{\infty}\frac{h(z')-h(z)}{(z'-z)^2}\,dz'.
\]

Its Fourier multiplier is `-|k|/2`. In particular,
`J[A cos(kz)] = -|k| A cos(kz)/2`.

The first-order in-plane crack-front weight-function relations in an infinite
isotropic solid are

\[
\begin{aligned}
\delta K_1&=K_1^0\mathcal J[h],\\
\delta K_2&=p_2 K_2^0\mathcal J[h]-b_2 K_3^0 h',\\
\delta K_3&=p_3 K_3^0\mathcal J[h]+b_3 K_2^0 h',
\end{aligned}
\]

where

\[
p_2=\frac{2-3\nu}{2-\nu},\qquad
p_3=\frac{2+\nu}{2-\nu},\qquad
b_2=\frac{2}{2-\nu},\qquad
b_3=\frac{2(1-\nu)}{2-\nu}.
\]

These are the planar crack's linear elastic response, including the change to
the local front frame. With the convention in the question,
`n = e_x - h' e_z + O(h'^2)` and `t = e_z + h' e_x + O(h'^2)`.
This fixes the signs of the slope terms. As a sign check, when `nu = 0` their
coefficients reduce to the local rotation of the two shear components:
`delta K2 = -K3^0 h'` and `delta K3 = K2^0 h'`, in addition to the line terms.
No term proportional to `h dK_alpha^0/da` occurs: the prescribed straight-front
intensities are independent of front position.

Write `r2 = K2^0/K1^0 = 0.6` and `r3 = K3^0/K1^0 = 0.8`, and set

\[
\begin{aligned}
D&=1+r_2^2+\frac{r_3^2}{1-\nu}
  =2.274285714285714,\\
H&=1+p_2r_2^2+\frac{p_3r_3^2}{1-\nu}
  =2.469915966386555,\\
L&=H/D=1.086018326928761.
\end{aligned}
\]

Linearizing Irwin's relation gives

\[
\frac{\delta G}{G^0}
=\frac{2}{D}\left(\frac{\delta K_1}{K_1^0}
 +r_2\frac{\delta K_2}{K_1^0}
 +\frac{r_3}{1-\nu}\frac{\delta K_3}{K_1^0}\right)
=2L\mathcal J[h]=-L|\partial_z|h.
\]

The slope terms cancel because `b3/(1-nu) = b2`.

## 2. Part A

For `h = A cos(2 pi z/lambda)`, `h' = -2 pi (A/lambda) sin(2 pi z/lambda)`.
Consequently,

\[
\begin{aligned}
C_1&=-\pi,& S_1&=0,\\
C_2&=-\pi p_2r_2,& S_2&=2\pi b_2r_3,\\
C_3&=-\pi p_3r_3,& S_3&=-2\pi b_3r_2,\\
C_G&=-2\pi L,& S_G&=0.
\end{aligned}
\]

Numerically the eight values, in this order, are

```
-3.141592653589793, 0,
-1.219677147864273, 5.913586171463141,
-3.400312048591305, -3.104632740018148,
-6.823654395086550, 0.
```

## 3. Part B

Both `J[f]` and `f'` have zero mean for a periodic front. Since `g` also has zero
mean, averaging equilibrium gives `G0 = Gamma0 + O(epsilon^2)`. Thus

\[
K_{10,\mathrm{scaled}}=D^{-1/2}=0.6630976492661631.
\]

The Poisson-kernel toughness has Fourier series

\[
g(s)=2\sum_{n=1}^{\infty}q^n\cos(2\pi ns).
\]

For `a = a_bar + epsilon lambda f(s)`, equilibrium at first order is
`-L |partial_s| f = g`. The zero-mean solution is therefore

\[
\begin{aligned}
f(s)&=-\frac{1}{\pi L}\sum_{n=1}^{\infty}
          \frac{q^n}{n}\cos(2\pi ns)\\
    &=\frac{\log d(s)}{2\pi L},\\
d(s)&=1-2q\cos(2\pi s)+q^2,\\
f'(s)&=\frac{2q\sin(2\pi s)}{L d(s)}.
\end{aligned}
\]

The logarithm has zero average for `q < 1`; hence the required mean convention
is already satisfied. The solution is unique up to the removed constant mode.
The front is furthest behind at the tough-band centers `s = 0 mod 1` and
furthest ahead at `s = 1/2 mod 1`. Its peak-to-peak displacement is

\[
\max f-\min f
=\frac{\log[(1+q)/(1-q)]}{\pi L}
=0.3220011521858956.
\]

The three SIF profiles used for the sixteen samples are

\[
\begin{aligned}
\kappa_1(s)&=\frac{g(s)}{2L},\\
\kappa_2(s)&=\frac{p_2r_2g(s)}{2L}-b_2r_3f'(s),\\
\kappa_3(s)&=\frac{p_3r_3g(s)}{2L}+b_3r_2f'(s).
\end{aligned}
\]

To find the continuous maximum of mode III, set `theta = 2 pi s` and write

\[
\kappa_3=\frac{a+b\sin\theta}{h-d\cos\theta}-c,
\quad c=\frac{p_3r_3}{2L},\quad a=c(1-q^2),\quad
b=\frac{2qb_3r_2}{L},\quad h=1+q^2,\quad d=2q.
\]

All stationary points satisfy
`b h cos(theta) - a d sin(theta) = b d`. There are two per period. Evaluating
both and selecting the larger gives

\[
\kappa_{3,\max}=1.231909965952339,
\qquad s_{\kappa_{3,\max}}=0.04092502626235328.
\]

The other stationary point is a minimum at approximately `s = -0.2259840864`.

## 4. Checks and limits

`solve.py` performs all checks before writing the result:

- Reconstructs `CG` and `SG` directly from the SIF coefficients and Irwin's
  relation, including cancellation of the energy sine coefficient.
- Checks the straight-front loading normalization, zero means, periodicity,
  and reflection symmetry of the front.
- Independently sums and solves 96 Fourier harmonics and compares every
  profile against the closed forms at 4096 positions. Since `q = 0.5`, omitted
  Fourier tails are below double-precision relevance.
- Evaluates the periodic real-space line integral by adaptive quadrature at
  all sixteen requested positions. Its kernel is obtained from
  `sum_m 1/(u+m)^2 = pi^2/sin(pi u)^2`; symmetrization removes the principal
  value singularity. This checks the normalization and sign of the line
  operator against its Fourier implementation.
- Reconstructs the pointwise first-order energy release rate from all three
  SIF profiles at 4096 positions and checks it against `g`.
- Checks the analytic peak-to-peak displacement against a grid containing both
  front extrema.
- Locates the largest mode III value on a 65536-point grid over a full period,
  refines it using scalar numerical optimization, and compares it with the
  analytic stationary-point calculation.
- Checks the requested JSON schema, finite numerical values, and sixteen
  entries per sampled profile.

All checks passed. The Fourier/profile and equilibrium residuals were below
`2e-14`; real-space integral discrepancies were below `2e-11`. The optimized
mode III maximum agreed with the analytic value within `1e-12`, and its location
within `1e-8` of a period.

These checks validate the first-order calculation and its numerical evaluation.
No independent three-dimensional elastic boundary-value computation was run,
and no finite-epsilon or higher-order corrections were evaluated.
