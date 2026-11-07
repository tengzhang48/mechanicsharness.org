# Energy Release for Periodic Film Cracks

## Results

For the given plane-stress, local-foundation model, define

\[
\bar E=\frac{E}{1-\nu^2},\qquad
\boxed{\ell=\sqrt{\frac{h\bar E}{k}}},\qquad
w=\frac H\ell,\qquad T=\tanh\frac w2.
\]

The energy releases per advancing tip, using the stipulated advancing area
normalization, are

\[
\boxed{G_A=\frac{\sigma_0^2\ell}{\bar E}\tanh\frac{H}{2\ell}},
\]
\[
\boxed{G_B=\frac{\sigma_0^2\ell}{\bar E}
\left[2\tanh\frac{H}{2\ell}-\tanh\frac H\ell\right]
=\frac{\sigma_0^2\ell}{\bar E}\frac{2T^3}{1+T^2}}.
\]

Both geometries are pure Mode I. In particular,

\[
K_{II}=0,\qquad
\boxed{\frac{K_{I,j}}{\sigma_0\sqrt\ell}
=\sqrt{(1-\nu^2)\gamma_j}},\qquad
\gamma_j=\frac{G_j\bar E}{\sigma_0^2\ell},\quad j=A,B.
\]

The positive root is selected for tensile crack opening. Mode III is absent
from this in-plane model. All requested quantities have been obtained;
there are no unresolved specifications or missing requested results.

## Strip Fields

Consider a translation-invariant relaxed strip of width W, centered at c.
Write y=x2-c, with -W/2<y<W/2. The bounded solution has u1=0. Equilibrium
for the normal displacement is

\[
h\bar E u_2''-ku_2=0.
\]

At both physical horizontal strip edges, traction freedom requires

\[
\sigma_{22}=0,\qquad\sigma_{12}=0,
\qquad u_2'(\pm W/2)=-\sigma_0/\bar E.
\]

No condition sets sigma11 to zero at these horizontal faces; it is tangential
stress, not traction on a horizontal edge. There is also no requirement of
displacement continuity across a crack. The exact solution is

\[
\boxed{u_2(y)=-\frac{\sigma_0\ell}{\bar E}
\frac{\sinh(y/\ell)}{\cosh(W/(2\ell))}},\qquad u_1=0,
\]
\[
\epsilon_{22}=-\frac{\sigma_0}{\bar E}
\frac{\cosh(y/\ell)}{\cosh(W/(2\ell))},\qquad
\epsilon_{11}=\epsilon_{12}=0,
\]
\[
\boxed{\sigma_{22}=\sigma_0\left[1-
\frac{\cosh(y/\ell)}{\cosh(W/(2\ell))}\right]},
\]
\[
\boxed{\sigma_{11}=\sigma_0\left[1-\nu
\frac{\cosh(y/\ell)}{\cosh(W/(2\ell))}\right]},
\qquad\sigma_{12}=0.
\]

The two edge displacements are opposite, and each equals in magnitude
(sigma0 ell/Ebar) tanh(W/(2ell)). The displacement is directed away from
each crack line. The residual edge stress is sigma11=sigma0(1-nu).

The length ell controls normal wake relaxation away from a crack face.
It is not a claim that every component of a two-dimensional tip field has
one screening length: the shear modulus is mu=E/[2(1+nu)], and a transverse
elastic disturbance has screening length sqrt(h mu/k). The latter is not
excited in these translation-invariant strip states. For u1, equilibrium
and shear-free edges give h mu u1''-k u1=0 and u1'=0 at both edges, whose
only solution is u1=0. Positive elastic and foundation quadratic energies
also give uniqueness of the strip equilibrium.

The far fields for the specified arrays are as follows.

| Case And Direction | Strip Width W | Strip Center c | Field |
| --- | ---: | --- | --- |
| A, x1 to +infinity | No physical strip edges | Not applicable | u=0, sigma=sigma0 I |
| A, x1 to -infinity | H | (m+1/2)H | Relaxed strip solution above |
| B, x1 to +infinity | 2H | (2m+1)H | Relaxed strip solution above |
| B, x1 to -infinity | H | (m+1/2)H | Relaxed strip solution above |

In A behind the tips, every edge x2=mH is a free crack face. In A ahead,
the same lines are intact, not traction-free edges: displacement and
traction are continuous, and reflection symmetry gives u2=0 and
sigma12=0 there. The full ahead field is the prestress state.

In B ahead, the edges x2=2mH and x2=(2m+2)H are free faces of old cracks;
the center x2=(2m+1)H is intact and has

\[
u_2=0,\quad \sigma_{12}=0,\quad
\sigma_{22}=\sigma_0[1-\operatorname{sech}(H/\ell)].
\]

In B behind, both old and new crack lines are free strip edges.
The specified periodicity maps corresponding strip interiors and face
traces; it does not remove the allowed jumps across crack lines.

## Energy Accounting

Let U(W) be the change in potential energy per unit x1 length of a relaxed
width-W strip. Substitution in the given functional gives

\[
U(W)=\int_{-W/2}^{W/2}
\left[h\sigma_0 u_2'+\frac{h\bar E}{2}(u_2')^2
+\frac k2u_2^2\right]dy.
\]

Multiplying h Ebar u2''-k u2=0 by u2 and integrating by parts yields

\[
\int_{-W/2}^{W/2}[h\bar E(u_2')^2+ku_2^2]dy
=h\bar E[u_2u_2']_{-W/2}^{W/2}
=-h\sigma_0[u_2]_{-W/2}^{W/2}.
\]

Consequently,

\[
\boxed{U(W)=\frac{h\sigma_0}{2}[u_2]_{-W/2}^{W/2}
=-\frac{h\sigma_0^2\ell}{\bar E}\tanh\frac W{2\ell}}.
\]

The first equality expresses half the prestress boundary work at elastic
equilibrium; it includes foundation energy, not just film strain energy.

To avoid infinities, truncate a periodic cell far ahead and far behind and
consider a small common tip advance da. Translating the steady equilibrium
transition region does not change its own integrated energy. The only net
change, in the limit of remote truncation boundaries, is replacing an ahead
slice by a behind slice. Hence

\[
G=\frac{U_{\rm ahead}-U_{\rm behind}}{Nh},
\]

where N is the number of advancing tips per transverse cell.

| Case | Period Width | Advancing Tips N | U Ahead | U Behind |
| --- | ---: | ---: | --- | --- |
| A | H | 1 | 0 | U(H) |
| B | 2H | 1 | U(2H) | 2U(H) |

Thus A releases -U(H)/h and B releases [U(2H)-2U(H)]/h, proving the
formulas stated above. An average energy density can alternatively be
multiplied by the corresponding period width; it must not be used directly
as energy per tip. The advancing area is h da, not 2h da. Existing cracks
in B are stationary and contribute no additional advancing tips. Their
changed opening and elastic energy are already included in U(2H) and 2U(H).

Reflection about each moving crack line preserves the geometry, loading
and foundation. Uniqueness gives u1 even and u2 odd under this reflection,
excluding antisymmetric Mode II. The stipulated local plane-stress
identity G=(KI^2+KII^2)/E then gives KI=sqrt(E G). The use of E here,
rather than Ebar, explains the factor 1-nu^2 in the normalized intensity.
This intensity is inferred from energy and symmetry, not from an independent
numerical fit to a crack-tip stress field.

## Critical Stress

For a constant positive toughness Gc and the prescribed advance pattern,

\[
\boxed{\sigma_{c,A}=\sqrt{\frac{\bar E G_c}
{\ell\tanh(H/(2\ell))}}},
\]
\[
\boxed{\sigma_{c,B}=\sqrt{\frac{\bar E G_c}
{\ell[2\tanh(H/(2\ell))-\tanh(H/\ell)]}}}.
\]

Equality gives the ideal quasistatic onset criterion; G>Gc supplies a
driving force for the specified advance but does not prescribe speed.
For all finite H/ell>0,

\[
0<\gamma_B<\gamma_A<1,\qquad
\frac{\gamma_B}{\gamma_A}=\frac{2T^2}{1+T^2}<1,
\qquad \sigma_{c,B}>\sigma_{c,A}.
\]

Thus a new crack halfway between existing cracks is more strongly shielded
than simultaneous cracks entering an initially uncracked film.

## Spacing Limits

Here h, E, nu and k are fixed, so ell is fixed.

For w=H/ell tending to zero,

\[
\gamma_A=\frac w2-\frac{w^3}{24}+\frac{w^5}{240}+O(w^7),
\qquad
\gamma_B=\frac{w^3}{4}-\frac{w^5}{8}
+\frac{17w^7}{320}+O(w^9).
\]

Therefore

\[
G_A\sim\frac{\sigma_0^2H}{2\bar E},\qquad
G_B\sim\frac{\sigma_0^2H^3}{4\bar E\ell^2},
\]
\[
\sigma_{c,A}\sim\sqrt{\frac{2\bar E G_c}{H}},\qquad
\sigma_{c,B}\sim2\ell\sqrt{\frac{\bar E G_c}{H^3}},
\]
\[
\frac{K_{I,A}}{\sigma_0\sqrt\ell}
\sim\sqrt{\frac{(1-\nu^2)H}{2\ell}},\qquad
\frac{K_{I,B}}{\sigma_0\sqrt\ell}
\sim\frac{\sqrt{1-\nu^2}}2(H/\ell)^{3/2}.
\]

In a narrow relaxed strip, sigma22 approaches zero, u2 approaches
-sigma0 y/Ebar, and sigma11 approaches sigma0(1-nu). There is substantial
normal relaxation but no far-field tangential strain. A releases the
approximately constant normal-relaxation energy density over a width H,
so its per-tip release is proportional to H. In B, the old width-2H strip
is already nearly relaxed; splitting it adds only an order-H^3 gain.
For example, the far-ahead normal stress on the future crack line is only
sigma0[1-sech(w)]~sigma0 w^2/2. Both critical stresses diverge as spacing
vanishes, and the divergence is stronger for B. Actual materials may
yield or invalidate small-strain assumptions before those divergences.

For H/ell tending to infinity,

\[
\gamma_A=1-2e^{-H/\ell}+O(e^{-2H/\ell}),\qquad
\gamma_B=1-4e^{-H/\ell}+O(e^{-2H/\ell}),
\]
\[
\boxed{G_A,G_B\longrightarrow\frac{\sigma_0^2\ell}{\bar E}},\qquad
\boxed{\sigma_{c,A},\sigma_{c,B}\longrightarrow
\sqrt{\frac{\bar E G_c}{\ell}}},
\qquad
\frac{K_{I,A}}{\sigma_0\sqrt\ell},
\frac{K_{I,B}}{\sigma_0\sqrt\ell}\longrightarrow\sqrt{1-\nu^2}.
\]

Wake relaxation is confined to distances of order ell from each face.
Widely separated cracks then interact negligibly, and both cases approach
the same isolated-crack release. For nu=0.3 the limiting intensity is
sqrt(0.91)=0.9539392014.

## Numerical Samples

The requested samples use nu=0.3. These are dimensionless values, not
choices of dimensional film or foundation properties.

| ell/H | GA Ebar/(sigma0^2 ell) | GB Ebar/(sigma0^2 ell) | KI,A/(sigma0 sqrt(ell)) | KI,B/(sigma0 sqrt(ell)) |
| ---: | ---: | ---: | ---: | ---: |
| 0.5 | 0.7615941560 | 0.5591607318 | 0.8324966558 | 0.7133276007 |
| 1 | 0.4621171573 | 0.1626401586 | 0.6484802334 | 0.3847109880 |
| 1.5 | 0.3215127375 | 0.0602425297 | 0.5409034952 | 0.2341382114 |
| 2 | 0.2449186624 | 0.0277201675 | 0.4720974293 | 0.1588249114 |

Dimensional units are ell in m, k in Pa/m, G in J/m^2, sigma_c in Pa,
and KI in Pa sqrt(m). U(W) is energy per x1 length, in J/m. All five
columns of the sample table are dimensionless. Full precision values are
saved in results.json.

## Executed Checks

The reproducible command is

```bash
timeout 120s python check.py
```

It uses only local Python libraries, logs PID and UTC time, and saves
actual outputs incrementally in CHECKS.txt and results.json. The executed
environment was Python 3.13.9, NumPy 2.4.3, SciPy 1.18.0, SymPy 1.14.0.
All eight checks passed with the following outputs and declared tolerances.

| Check | Actual Output | Tolerance | Errors Detected |
| --- | --- | --- | --- |
| Symbolic equilibrium and both normal face tractions; positive B identity | All four residuals exactly 0 | Exact zero | Wrong length, displacement sign, face loading, or B subtraction |
| Symbolic small-spacing series | A: w/2-w^3/24+w^5/240-17w^7/40320; B: w^3/4-w^5/8+17w^7/320 | Exact displayed coefficients | Wrong dense-array powers or coefficients |
| Direct integration of all three energy terms over six strip widths | Maximum absolute normalized-energy error 1.1102230246251565e-16 | 1e-11 | Missing foundation term, wrong prestress work, factor of two |
| Numerical boundary-value solve with zero initial guess | Maximum field error 1.5587531265737198e-13; face residual 0 | Field 2e-9; face 1e-12 | Wrong strip centering, sign or edge condition |
| Independently assembled variational finite elements | Final max relative error A 1.7824124855032374e-8; B 1.8531206615507827e-7 | Both below 2e-6 | Wrong analytic energy or old/new strip accounting |
| Small and large spacing limits | At w=0.001, A/(w/2)=0.999999916666675, B/(w^3/4)=0.9999995000002124; at w=40, both gamma=1 | Small relative 1e-6; large absolute 1e-12 | Wrong asymptotic powers or isolated normalization |
| Plane-stress KI conversion | Maximum gamma reconstruction error 5.551115123125783e-17 | 1e-14 | E versus Ebar normalization error, not independent KI validation |
| Shielding and critical-stress ordering | All 0<B<A<1; sigma_c,B/sigma_c,A=1.1670607655, 1.6856296119, 2.3101888922, 2.9724394322 | Strict inequalities | Interchanged cases or release sign errors |

For the variational route, dimensionless displacement v is represented
with continuous linear elements. The assembled quadratic functional is
integral(v'+v'^2/2+v^2/2) dz, with natural boundary loads obtained directly
from its first variation. The computed minimizer is not filled from the
analytic hyperbolic solution. Energy is then integrated element by element.

| Elements Per Strip | Maximum Relative Error A | Maximum Relative Error B |
| ---: | ---: | ---: |
| 128 | 4.5627680228e-6 | 4.7436285673e-5 |
| 256 | 1.1407311496e-6 | 1.1859749658e-5 |
| 512 | 2.8518523381e-7 | 2.9649798079e-6 |
| 1024 | 7.1296461373e-8 | 7.4124760166e-7 |
| 2048 | 1.7824124855e-8 | 1.8531206616e-7 |

Successive B-error reduction factors were 3.9997712465, 3.9999428078,
3.9999857016, and 3.9999964224, within the declared range 3.5 to 4.5.
This is the expected second-order energy convergence.

The numerical BVP and variational minimization use distinct numerical
constructions but share the stated elastic/foundation law and the
translation-invariant strip assumption. Quadrature independently checks
integration and bookkeeping, but uses the analytic strip displacement.
The symbolic formulas, their series and the KI algebra are consistency
checks, not independent evidence for the physical model. No identical
formula evaluated twice is presented as an independent solution. No
two-dimensional crack-tip stress extraction or experimental check was run.

## Model Scope

Supported conclusions within this model are the exact steady per-tip
releases for the imposed equal-advance patterns, pure Mode I by symmetry,
the ideal toughness thresholds, exponential screening at wide spacing,
and stronger shielding of infill cracks in B. Old cracks can change their
opening while remaining fixed, and a relaxed strip retains tangential
stress.

The model does not determine crack nucleation, selected spacing, arbitrary
unequal-tip stability, growth kinetics, or evolution of a general crack
network. It does not establish validity of a local foundation for a real
three-dimensional substrate, nor include substrate nonlocality, plasticity,
interfacial debonding, bending, out-of-plane fracture, or a finite process
zone. Applying the thresholds to materials requires dimensional h, E, k,
sigma0, toughness data and validation of linear elasticity, plane stress,
bond maintenance and the local-foundation approximation.

No scientific sources, other runs, or reference solutions were consulted.
No necessary specification was unclear, so QUESTIONS.md is not required.
Missing requested work: none. Independent 2D tip extraction and empirical
validation remain explicitly outside the executed checks.
