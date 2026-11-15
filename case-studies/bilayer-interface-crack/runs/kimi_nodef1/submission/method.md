# Method

## Exact energy and amplitude relations

I used the plane-strain modulus \(\bar E_i=E_i/(1-\nu_i^2)\). Because the material far ahead of the tip is stress-free, the steady-state energy release rate equals the strain energy per unit length in the two independent beams far behind it:

\[
G=\frac{1}{2\bar E_1}\left(\frac{P^2}{h}+\frac{12M^2}{h^3}\right)
+\frac{1}{2\bar E_2}\left[\frac{P^2}{H}+
\frac{12\{M+P(h+H)/2\}^2}{H^3}\right].
\]

Here \(H=2h\), \(\bar E_1=(81/8)E_2\), and \(\bar E_2=(9/8)E_2\).
Consequently,

\[
\frac{GE_2h}{P^2}\bigg|_{M=0}=\frac{287}{162},\qquad
\frac{GE_2h^3}{M^2}\bigg|_{P=0}=\frac{34}{27}.
\]

For the complex intensity convention used below,

\[
T(r)\sim \frac{K}{\sqrt{2\pi r}}r^{i\epsilon},
\qquad A=|K|.
\]

With plane strain \(\kappa_i=3-4\nu_i\), the Dundurs parameters are
\(\alpha=4/5\), \(\beta=1/5\), and

\[
\epsilon=\frac{1}{2\pi}\log\frac{1-\beta}{1+\beta}
=-0.06453177620670408.
\]

Also
\[
D\equiv\frac{\kappa_1+1}{\mu_1}+\frac{\kappa_2+1}{\mu_2}
=8\left(\frac1{\bar E_1}+\frac1{\bar E_2}\right)
=\frac{640}{81E_2}.
\]
The Irwin closure integral gives
\[
G=\frac{D}{16\cosh^2(\pi\epsilon)}A^2
=\frac{64}{135}\frac{A^2}{E_2},
\]
where \(\operatorname{sech}^2(\pi\epsilon)=1-\beta^2=24/25\).
Thus the two nondimensional amplitudes are exactly
\(\sqrt{1435/384}\) and \(\sqrt{85/32}\), respectively.
The decade phase change is \(\epsilon\log 10\), converted to degrees.

## Numerical determination of the complex intensity phases

The energy relation fixes \(|K|\), but not its phase. I therefore solved the
plane-strain problem with quadratic finite elements. The implementation is in
`fem_solve.py` and `fem_quadratic.py`; the saved displacement data are in the
`q9_*.json` files.

I set \(h=1\), \(H=2\), \(E_2=1\), and \(E_1=9\). The computational domain was
\([-L,L]\times[-H,h]\), with \(L=20h\) for the final fit. The top and bottom
outer surfaces and the right end were traction-free. The specified beam
stresses were applied as consistent nodal tractions on the left ends. The
negative interface was split into two traction-free crack faces, while the
positive interface had common displacements. Three constraints only removed
rigid-body motion; the self-equilibrated applied loads make their reactions
zero.

A polar tensor-product mesh was geometrically graded toward the crack tip.
Nine-node quadrilaterals and six-node center triangles were used, with 3 by 3
and three-point quadrature, respectively. The final mesh (`q9_120_L20_n128`)
had 124,781 nodes, a minimum radial scale of \(10^{-5}h\), and 128 nominal
angular subdivisions per layer (with additional rays through the rectangular
corners). Cases \(P=1,M=0\) and \(P=0,M=1\) were solved separately and combined
by superposition.

For a distance \(s\) behind the tip, let the upper-minus-lower displacement
jump be \(\delta=\Delta u_2+i\Delta u_1\). In the same intensity convention,
\[
\delta(s)=rac{D}{2(1+2i\epsilon)\cosh(\pi\epsilon)}
K\sqrt{\frac{s}{2\pi}}s^{i\epsilon}.
\]
I evaluated the left side at the crack-face nodes and inverted this equation
for a raw \(K(s)\). It was extrapolated by the complex least-squares model
\[
K(s)=K_0+a_1s^{1/2-i\epsilon}+a_2s+O(s^{3/2-i\epsilon})
\]
over \(0.001\le s/h\le0.02\). The powers account for the regular
\(O(s)\) opening contribution and the next \(O(s^{3/2+i\epsilon})\)
singularity-family contribution. Varying the fit interval and adding the
latter term explicitly changed the reported phases by only a few hundredths
of a degree.

Because lengths were measured in units of \(h\), the oscillatory factor is one
at \(r=h\); hence `phase_at_h_deg` is simply \(\arg K_0\) in degrees. For the
combined load \(M=Ph\), I added the two unit-load intensity vectors. To avoid
retaining the small finite-element amplitude error, each direction from the
FE fit was assigned its exact energy-derived amplitude before this addition.
`produce_results.py` contains this final calculation and writes `results.json`.

## Checks

- The final FE/exact amplitude ratios were **0.9990032** for case A and
  **0.9971887** for case B before replacing their amplitudes by the exact
  energy values. This independently checks the opening-jump prefactor, load
  normalization, and energy balance.
- Increasing radial refinement from 150 to 200 intervals at otherwise
  identical settings changed the fitted phases from 52.29515, 4.34261 degrees
  to 52.29486, 4.34201 degrees. The final 128-angle refinement changed case A
  and case B by about 0.003 and 0.006 degrees from the 96-angle calculation.
- Increasing \(L/h\) from 10 to 20 changed the phases by 0.027 degree (A) and
  0.059 degree (B), using the 96-angle mesh. Fit-window changes contributed
  less than about 0.04 degree in the final calculations.
- The final scaled sparse-solve relative residuals on the unconstrained
  equations were \(1.3\times10^{-10}\) (A) and \(4.6\times10^{-10}\) (B).
  Applied end-load force and moment resultants were zero to roundoff.
- A separate equal-layer pure-moment DCB calculation returned an amplitude
  3.44219, compared with the homogeneous analytical value
  \(\sqrt{12}=3.46410\) (0.63% difference); see `check_homogeneous.py`.

The energy rates, amplitudes, and decade phase change are analytical. The
three absolute phases remain numerical results; their mesh, fit-window, and
finite-domain tests are checks of convergence rather than independent exact
phase identities.
