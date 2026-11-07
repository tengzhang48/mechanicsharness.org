# Steady viscoelastic strip crack: solution and verification

## Result and reproduction

There is one steady state at every speed

\[
\boxed{0<V<V_*},\qquad \boxed{V_*=13.4235715733}.
\]

The endpoint is a supremum, **not** a state with a finite cohesive zone.
As it is approached, \(\epsilon\to\infty\), \(a\to1\), and
\(R\to5\). The physical endpoint speed is \(v_*=V_*h/(\pi t_c)\).
All three requested speeds are admissible.

Run from the task directory:

```sh
python solve.py case.json results.json
```

`results.json` contains all 45 requested field rows, the scalar results,
spatial-integral data, convergence and tail checks, and actual verification
residuals. `run.log` records the submitted run. Only Python, NumPy and SciPy
are required. The code saves partial JSON before each case's verification.
Optional evidence assertions can be reproduced with
`python verify.py case.json results.json`; `verification.log` contains the
successful output. These assertions check consistency, not an external
reference solution.

The numerical tables in section 6 summarize the submitted run. Numerical
bulk loss is always calculated by spatial integration, never by setting it
equal to \(R-1\).

## 1. Material reduction and stress potential

Define the dimensionless displacement

\[
H=\frac{\pi G_\infty w}{\tau_ch},\qquad
(H_X,H_Y)=(E_1,E_2).
\]

A material point travels towards decreasing \(X\), so
\(\dot f=-(V/t_c)f_X\). The stated creep law, including the fully relaxed
initial history, gives

\[
E_i(X,Y)=\lambda T_i(X,Y)+(1-\lambda)\mathcal Q_VT_i(X,Y),\qquad
\mathcal Q_V f=\frac1V\int_X^\infty e^{-(s-X)/V}f(s,Y)\,ds.
\tag{1}
\]

In particular, a constant far-ahead stress gives \(E_i=T_i\), not
\(E_i=\lambda T_i\). Equivalently,

\[
E_i-VE_{i,X}=T_i-\lambda VT_{i,X}.
\tag{2}
\]

Strain compatibility and (2) imply
\((1-\lambda V\partial_X)(T_{2,X}-T_{1,Y})=0\).
The possible homogeneous term grows as \(e^{X/(\lambda V)}\) and is excluded
by the relaxed bounded far-ahead condition. Thus stress is curl-free.
Together with equilibrium this allows a harmonic potential \(U\):

\[
(T_1,T_2)=(U_X,U_Y),\quad \nabla^2U=0,\quad
H=\lambda U+(1-\lambda)\mathcal Q_VU.
\tag{3}
\]

On the top edge \(H_X=0\); equation (2) and the bounded far-ahead history
give \(T_1=0\) there. Similarly, \(T_1=0\) on the intact interface. Fixing
the potential constant with the far-ahead state gives

\[
U(X,\pi)=\pi a,\quad U(X\geq0,0)=0,\quad
U_Y(X,0)=\begin{cases}1&-\epsilon<X<0,\\0&X<-\epsilon.\end{cases}
\tag{4}
\]

The stress problem is therefore elastic-looking, but the displacement is
the creep transform (3), not \(U\) itself.

## 2. Explicit stress and displacement fields

Let

\[
z=X+iY,\quad q=e^z,\quad b=e^{-\epsilon},\quad
L=\sqrt{e^\epsilon-1},\quad B(q)=\frac{\sqrt q}{\sqrt{q-1}}.
\]

The square roots continue from the upper half \(q\)-plane; hence
\(B=-i\sqrt{q/(1-q)}\) on the cracked/cohesive bottom edge. Write
\(S=T_2+iT_1\). A Cauchy representation satisfying the prescribed real
traction on that edge is

\[
S(q)=B(q)\left[a+\frac1\pi\int_b^1
\frac{\sqrt{(1-t)/t}}{t-q}\,dt\right].
\tag{5}
\]

On \(b<q<1\), the imaginary part of the bracket supplies
\(\Re S=1\); on \(0<q<b\), \(\Re S=0\). On the top and intact boundaries
\(S\) is real. At \(q\to\infty\), \(S\to a\).
Removal of the \((q-1)^{-1/2}\) front singularity requires the bracket
in (5) to vanish at \(q=1\):

\[
\boxed{a(\epsilon)=\frac1\pi\int_b^1\frac{dt}{\sqrt{t(1-t)}}
=\frac2\pi\arccos(e^{-\epsilon/2})=\frac2\pi\arctan L.}
\tag{6}
\]

Substituting \(t=\cos^2\theta\) in (5) gives the elementary field

\[
\boxed{S(z)=\frac2\pi\arctan\left(\frac{L}{\sqrt{1-e^{-z}}}\right),
\quad T_1=\Im S,\quad T_2=\Re S.}
\tag{7}
\]

The arctangent uses the continuation with \(\Re S=1\) on the cohesive
edge, not the opposite logarithm branch. The implementation evaluates
equivalent expressions on either side of \(X=0\) to avoid exponential
overflow. It uses explicit real-edge formulas on \(Y=0\).

A convenient potential, including its additive constant, is

\[
\boxed{U(X,Y)=\pi a-\int_Y^\pi\Re S(X+i\eta)\,d\eta.}
\tag{8}
\]

Analyticity of \(S\) and \(T_1(X,\pi)=0\) give \(U_X=T_1\).
Equations (1), (3), (7), and (8) specify **all** stress, strain and
displacement fields. Physical values follow from
\(\tau_i=\tau_cT_i\), \(\gamma_i=(\tau_c/G_\infty)E_i\), and
\(w=(\tau_ch/(\pi G_\infty))H\). In the lower half, \(T_1,E_1,w\)
are odd in \(y\), while \(T_2,E_2\) are even.

The front is bounded: \(S\to1\). Ahead,
\(T_2=1-2\sqrt X/(\pi L)+O(X^{3/2})\); on the cohesive side,
\(T_1=-2\sqrt{-X}/(\pi L)+O((-X)^{3/2})\).
There is an allowed logarithmic stress singularity at the *rear* edge
\(X=-\epsilon,Y=0\), owing to the Dugdale traction drop. No point sample
is at that singularity, and its squared bulk loss is spatially integrable.

## 3. Opening equation, uniqueness, and speed range

Use \(r=X+\epsilon\) in the cohesive zone, so that its rear is \(r=0\)
and front is \(r=\epsilon\). Define

\[
t_\epsilon(r)=\frac2\pi\operatorname{atanh}
\sqrt{\frac{e^{-r}-e^{-\epsilon}}{1-e^{-\epsilon}}}
=-T_1(-\epsilon+r,0),\qquad
U(-\epsilon+r,0)=\int_r^\epsilon t_\epsilon(s)\,ds.
\tag{9}
\]

Direct integration yields
\(U(-\epsilon,0)=\pi a(\epsilon)^2/2\). This identity can also be
proved by differentiating both sides in \(\epsilon\) and using their
common zero limit at \(\epsilon=0\).

Applying (3) to the bottom displacement and interchanging the two
integrations gives the required scalar equation:

\[
\boxed{\mathcal H(\epsilon,V):=
\int_0^\epsilon[1-(1-\lambda)e^{-r/V}]t_\epsilon(r)\,dr
=\frac{\pi\bar\delta_c}{2}.}
\tag{10}
\]

For numerical stability at small \(V\), the code instead evaluates
\[
\mathcal H=\frac{\pi a^2}{2}-(1-\lambda)
\int_0^\epsilon e^{-r/V}t_\epsilon(r)\,dr.
\tag{11}
\]
The exponentially small omitted kernel beyond \(r=45V\) is negligible.
Logarithms are evaluated using `expm1` and `log1p` near the rear edge.

For fixed \(V>0\), \(t_\epsilon(r)\) and its domain increase with
\(\epsilon\), so \(\mathcal H\) strictly increases from zero.
For fixed \(\epsilon\), \(\mathcal H\) strictly decreases with \(V\).
Consequently a finite steady solution, if it exists, is unique, and
\(\epsilon,a,R\) increase with speed. The opening is nonnegative and
increases towards the rear because the bottom \(E_1\) is nonpositive.
Also, the intact traction decreases from 1 to \(a<1\), respecting the
pre-separation Dugdale condition.

As \(\epsilon\to\infty\),
\[
t_\infty(r)=\frac2\pi\operatorname{atanh}(e^{-r/2}),\qquad
U_\infty(r)=\frac4\pi\sum_{k=1,3,5,\ldots}\frac{e^{-kr/2}}{k^2}.
\]
Thus
\[
\boxed{\mathcal H_\infty(V)=\lambda\frac\pi2+
(1-\lambda)\frac4\pi\sum_{k=1,3,5,\ldots}
\frac1{k^2(1+kV/2)}.}
\tag{12}
\]

This decreases strictly from \(\pi/2\) to \(\lambda\pi/2\). For the
present \(\lambda<\bar\delta_c<1\), the finite endpoint solves
\(\mathcal H_\infty(V_*)=\pi\bar\delta_c/2\), giving the value above.
A finite \(\epsilon\) always has \(\mathcal H<\mathcal H_\infty\),
so \(V=V_*\) and all larger speeds are inadmissible.

More precisely, with
\(c_*=-\mathcal H_\infty'(V_*)=0.01024093897\),
\[
\mathcal H_\infty(V)-\mathcal H(\epsilon,V)
=2e^{-\epsilon/2}+o(e^{-\epsilon/2}),
\]
\[
\epsilon=-2\log(V_*-V)+2\log(2/c_*)+o(1),\qquad
1-a=\frac{c_*}{\pi}(V_*-V)+o(V_*-V).
\tag{13}
\]
The rear recedes to \(X=-\infty\) in the front-based coordinates.
At fixed interior \((X,Y)\), the endpoint field tends to
\((T_1,T_2,E_1,E_2)=(0,1,0,1)\). The nontrivial rear-attached limiting
opening profile is the \(U_\infty\) profile in (12), with the same
creep operation as before.

The slow-speed limit is
\[
\boxed{a_0=\sqrt{\bar\delta_c}=0.447213595500,\quad
\epsilon_0=-2\log\cos\left(\frac{\pi\sqrt{\bar\delta_c}}2\right)
=0.540362288055,\quad R_0=1.}
\tag{14}
\]
Away from boundary singularities \(E_i\to T_i\), and bulk dissipation
tends to zero. For example its leading global loss is
\[
R-1=\frac{2(1-\lambda)V}{\pi^2\bar\delta_c}
\left[\log\frac{4(1-e^{-\epsilon_0})}{V}+\gamma_{\rm Euler}\right]+o(V).
\]

## 4. Energy release and the local bulk-loss integrand

The Maxwell branch is relaxed ahead. The two half-strips therefore have
stored energy per unit length \(G_\infty\Delta^2/h\); behind both strain
and stress vanish. The definition of \(\Gamma\) in the question gives

\[
\boxed{\Gamma=\frac{G_\infty\Delta^2}{h},\qquad
R=\frac{a^2}{\bar\delta_c}.}
\tag{15}
\]

For each component set \(q_i=\tau_i-G_\infty\gamma_i\).
Equation (2) in material time implies
\[
\dot q_i+\frac{q_i}{\lambda t_c}
=(G_0-G_\infty)\dot\gamma_i.
\]
Hence the **relaxation** time is \(\lambda t_c\), not the specified
creep time \(t_c\), and the Maxwell viscosity is
\(\eta=(G_0-G_\infty)\lambda t_c=G_\infty(1-\lambda)t_c\).
Differentiating the supplied free energy gives exactly

\[
\boxed{\phi_D=\sum_i\frac{q_i^2}{\eta}
=\frac{\tau_c^2}{G_\infty(1-\lambda)t_c}
\big[(T_1-E_1)^2+(T_2-E_2)^2\big]\geq0.}
\tag{16}
\]

For a steady translating field the dissipation per unit advance is
\(v^{-1}\int\phi_D\,dx\,dy\). Including both halves and the coordinate
Jacobian gives the direct spatial integral used in the code:

\[
\boxed{\frac{\Gamma_D}{\Gamma_0}
=\frac{2}{\pi V(1-\lambda)\bar\delta_c}
\int_{-\infty}^{\infty}\int_0^\pi
\big[(T_1-E_1)^2+(T_2-E_2)^2\big]dY\,dX.}
\tag{17}
\]

The interface work is \(\tau_c\delta_c=\Gamma_0\). The global check is
therefore \(\Gamma_D/\Gamma_0+1=R\). Equation (17), not this check,
generates the reported loss. At the endpoint the energy balance predicts
\(\Gamma_D/\Gamma_0\to4\); that *limiting prediction* is distinguished
from the three independently integrated finite-speed values in the JSON.

## 5. Numerical method and evidence definitions

Scalar equations use adaptive quadrature and bracketed roots. Samples use
the analytic stress (7) and adaptive integration of (1). A second sample
calculation integrates \(e^{-t}T(X+Vt,Y)\) directly, without subtracting
the far-ahead state. Its agreement checks both history orientation and
far-state handling.

For the spatial bulk integral:

* Use core \(-\epsilon-32\leq X\leq32\), \(0\leq Y\leq\pi\).
* Use Gauss-Legendre quadrature after \(Y=\pi t^2\), resolving the lower
  edge. No singular region is removed.
* Add geometric \(X\)-nodes on each side of \(X=-\epsilon,0\), down to
  distance \(10^{-8}\), to resolve rear logarithms and front square roots.
* Approximate the analytic stress by cubic Hermite polynomials using its
  exact derivative. Integrate each polynomial against the exponential
  kernel exactly, propagating the history backwards in \(X\). Compare
  this spatial-grid history against independent adaptive quadrature at
  15 interior points per case.
* Integrate the nonnegative \(Y\)-summed squared density in \(X\) by the
  composite trapezoidal rule. For \(X<X_{\min}\), stress is already
  exponentially small. Retaining the homogeneous history
  \(Q(X)=Q(X_{\min})e^{(X-X_{\min})/V}\) gives the explicit raw tail
  \((1-\lambda)^2(V/2)\int_0^\pi |Q(X_{\min},Y)|^2dY\).
  The far-right constant is \(T_2=E_2=a\), contributing zero loss.
* Repeat with core lengths 24 and 40 in place of 32. This tests both the
  neglected stress correction in the left tail and the right truncation.

Four combined refinements use \((\text{maximum }\Delta X,N_Y,N_g)\)
equal to \((.02,48,100),(.01,80,160),(.005,128,240),(.0025,192,360)\),
where \(N_g\) is the geometric node count on each side of an interface end.
Additional checks refine \(X\), \(Y\), and geometric nodes separately
from the third grid. The final submitted loss uses the fourth grid.

Verification residuals use normalized components directly:

| Check | Residual or comparison | Errors it can detect |
|---|---|---|
| Equilibrium | \(T_{1,X}+T_{2,Y}\), centered finite differences | Wrong complex-component order, signs, coordinate mapping, nonharmonic field |
| Strain compatibility | \(E_{2,X}-E_{1,Y}\), centered finite differences | Inconsistent component histories, nonintegrable displacement |
| Material law | \(E-VE_X-T+\lambda VT_X\), centered finite differences | Reversed history, wrong modulus ratio, creep vs relaxation-time error |
| Local energy | \(\tau\dot\gamma-\dot\psi\) vs (16), finite-difference derivatives | Wrong material-time sign, free-energy coefficient, viscosity or component sum |
| Cohesive/crack traction | Bottom \(T_2=1/0\) at several points | Incorrect loading representation, branch or zone location |
| Intact displacement | Reconstruct \(U(X,0)\) vertically for three \(X>0\) | Incorrect potential constant or mixed-boundary solution; (3) then gives \(H=0\) |
| Top edge | \(E_1=0\) at three points; numerical history of \(U=\pi a\) | Drift of outer displacement or incorrect far-history normalization |
| Critical separation | Vertical reconstruction (8), followed by displacement creep, vs \(\bar\delta_c\) | Wrong factor of two or \(\pi\), stress mistaken for strain, or erroneous opening kernel |
| Front boundedness | Stress sequence approaching \(X=0\) from both sides | Residual square-root stress-intensity term |
| Far fields | Stress/strain at specified distant points | Missing relaxed offset, nonzero homogeneous constitutive mode, tail errors |
| Bulk energy balance | \((D+1-R)/\max(1,R)\) | Wrong Jacobian, omitted lower half, insufficient integration or omitted strain tail |

Derivative residuals have the dimensions removed by differentiating in
\(X,Y\): equilibrium is scaled by \(\tau_c\pi/h\), compatibility by
\(\tau_c\pi/(G_\infty h)\), and the material residual is scaled by
\(\tau_c/G_\infty\). The local energy residual is scaled by
\(\tau_c^2/(G_\infty t_c)\). Displacement checks use \(H,U\) as defined
above. These are absolute residuals, with no division by a small field.
Finite-difference steps \(.002,.001,.0005\) give approximately second-order
decay, rather than merely checking symbolic identities.

The independent vertical critical-opening calculation uses neither the
bottom-slope formula nor the integration-by-parts scalar equation (11).
Some boundary conditions are imposed exactly by the analytic construction;
their zero errors are not claimed to measure discretization accuracy.

## 6. Submitted numerical results

The scalar and convergence tables below are completed from the final run.
Every field row is explicitly stored under `cases[*].fields` in
`results.json`, in input coordinate order, with exactly the normalization
requested in `case.json`.

### Steady states and independently integrated bulk loss

| \(V\) | \(\epsilon\) | \(a\) | \(R\) | \(\Gamma_D/\Gamma_0\) | \(D+1-R\) |
|---:|---:|---:|---:|---:|---:|
| 0.001 | 0.544705555093 | 0.448842489170 | 1.007297900421 | 0.007299119833 | \(1.2194\times10^{-6}\) |
| 0.1 | 0.740780931670 | 0.514812186740 | 1.325157938082 | 0.325162283343 | \(4.3453\times10^{-6}\) |
| 10 | 6.333305684440 | 0.973162076321 | 4.735222133951 | 3.735222386309 | \(2.5236\times10^{-7}\) |

The relative energy discrepancies \(|D+1-R|/\max(1,R)\) are respectively
\(1.211\times10^{-6},3.279\times10^{-6},5.329\times10^{-8}\).
The displayed bulk-loss digits identify the quadrature output; they are
not a claim that every displayed digit is converged.

### Bulk-loss refinement

| Maximum \(\Delta X\), \(N_Y\), \(N_g\) | \(D\), \(V=.001\) | \(D\), \(V=.1\) | \(D\), \(V=10\) |
|---|---:|---:|---:|
| .02, 48, 100 | 0.007317669252 | 0.325295421383 | 3.735233542278 |
| .01, 80, 160 | 0.007305006578 | 0.325202481622 | 3.735225392151 |
| .005, 128, 240 | 0.007300883801 | 0.325172722056 | 3.735223010838 |
| .0025, 192, 360 | 0.007299119833 | 0.325162283343 | 3.735222386309 |

Last-refinement absolute changes are \(1.764\times10^{-6}\),
\(1.044\times10^{-5}\), and \(6.245\times10^{-7}\).
Separately refining the third grid gives a maximum change over all speeds
of \(9.409\times10^{-6}\) for \(X\), \(1.71\times10^{-12}\) for \(Y\),
and \(4.623\times10^{-6}\) for corner-node density. The full separate
checks, rather than only these maxima, are in the JSON.

The raw integrals below are of \(|T-E|^2\,dX\,dY\), before multiplying
by the factor in (17):

| \(V\) | Core \(X\) interval | Core nodes \(N_X\) | Raw core integral | Raw left-tail integral | Left-tail contribution to \(D\) |
|---:|---|---:|---:|---:|---:|
| .001 | [-32.5447055551, 32] | 27261 | \(2.0637775121\times10^{-6}\) | \(2.7454431015\times10^{-18}\) | \(9.7100\times10^{-15}\) |
| .1 | [-32.7407809317, 32] | 27340 | 0.00919374696517 | \(3.7862634506\times10^{-16}\) | \(1.3391\times10^{-14}\) |
| 10 | [-38.3333056844, 32] | 29577 | 10.5383440467 | 0.0227484408370 | 0.00804561512631 |

All three use the full \(0\leq Y\leq\pi\) range. Maximum stress magnitudes
on the left core boundary are \(4.64,5.18,7.16\) times \(10^{-8}\).
Maximum right-boundary stress deviations from \((0,a)\) are below
\(2.2\times10^{-15}\).

| Core length parameter \(B\): \([-\epsilon-B,B]\) | \(D\), \(V=.001\) | \(D\), \(V=.1\) | \(D\), \(V=10\) |
|---|---:|---:|---:|
| 24 | 0.007299119854 | 0.325162282303 | 3.735222799474 |
| 32 | 0.007299119833 | 0.325162283343 | 3.735222386309 |
| 40 | 0.007299119839 | 0.325162283985 | 3.735222383057 |

The maximum absolute tail/bounds changes relative to \(B=32\) are
\(2.095\times10^{-11},1.040\times10^{-9},4.132\times10^{-7}\).
At \(V=10\), changing \(B\) from 24 to 40 changes the explicit tail
contribution from 0.0398495919 to 0.00162438268 while leaving the total
stable. Thus the long strain wake is accounted for rather than discarded.

### Field and boundary verification

Here are maximum absolute residuals at the smallest finite-difference step
\(d=.0005\); larger-step data are saved in `results.json`.

| \(V\) | Equilibrium | Strain compatibility | Material law | Local power/free-energy identity | Critical \(\bar\delta_c\) error |
|---:|---:|---:|---:|---:|---:|
| .001 | \(3.705\times10^{-7}\) | \(1.681\times10^{-7}\) | \(1.667\times10^{-10}\) | \(5.658\times10^{-14}\) | \(1.50\times10^{-15}\) |
| .1 | \(3.732\times10^{-7}\) | \(3.081\times10^{-7}\) | \(1.464\times10^{-8}\) | \(8.660\times10^{-10}\) | \(3.49\times10^{-13}\) |
| 10 | \(3.691\times10^{-7}\) | \(2.113\times10^{-8}\) | \(3.473\times10^{-8}\) | \(3.560\times10^{-8}\) | \(1.56\times10^{-13}\) |

These residuals decrease by about four when the step is halved. The
spatial-grid creep calculation differs from independent adaptive
quadrature by at most \(3.85\times10^{-13}\) in a normalized component.
The alternative direct history kernel agrees with every requested strain
sample to better than \(2.2\times10^{-14}\).

Bottom cohesive and separated traction checks, and top \(E_1\), return
zero. Vertically reconstructed intact \(U\) is below
\(4.5\times10^{-16}\). The top numerical displacement history equals
\(\pi a\) to floating-point precision. Far-ahead normalized stress and
strain errors at \(X=32,Y=.4\pi\) are below \(2\times10^{-15}\).
Far-behind stress at \(X=-\epsilon-40,Y=.4\pi\) is below
\(1.1\times10^{-9}\); strain is checked farther back if necessary, at
\(X=-\epsilon-\max(40,40V)\), and is below \(8.1\times10^{-10}\).

For example, the requested front-near samples \((X,y/h)=(0,.1)\) are:

| \(V\) | \(T_1\) | \(T_2\) | \(E_1\) | \(E_2\) |
|---:|---:|---:|---:|---:|
| .001 | -0.218750363180 | 0.664451252249 | -0.218235438288 | 0.664361887719 |
| .1 | -0.192154239429 | 0.729800667772 | -0.157551089302 | 0.713539673712 |
| 10 | -0.009752727719 | 0.988576784492 | -0.001484025839 | 0.975883177790 |

The other 42 rows are in `results.json`, including all top-edge samples.

### Endpoint checks

The independent odd-series evaluation of (12), using 1024, 4096, and
16384 terms at the reported \(V_*\), has absolute opening residuals
\(1.018\times10^{-8},6.360\times10^{-10},3.975\times10^{-11}\).
Evaluating the limiting opening at \(V_*-.001\) and \(V_*+.001\) gives
signed differences from the target of \(+1.02416\times10^{-5}\) and
\(-1.02403\times10^{-5}\). At \(V=V_*\), finite lengths
\(\epsilon=8,16,24,32,40\) all fall short of the required opening, with
shortfalls \(1.8933\times10^{-2},4.9218\times10^{-4},
1.0484\times10^{-5},2.0686\times10^{-7},3.9386\times10^{-9}\).
This independently supports the finite supremum and infinite-length limit.

## 7. Scope, resources, and limitations

No physical input or requested finite-speed quantity is unresolved.
The endpoint's infinite length is represented as `null` in the auxiliary
`endpoint.epsilon_limit` field because JSON must contain finite numbers;
it is analytically infinite, not an unknown numerical result.

Errors reported are empirical convergence and residual estimates, not
interval-arithmetic certification. The ideal Dugdale rear-edge singularity
is retained, not regularized. The prescribed quasistatic linear model is
used without inertia or an additional fracture law.

Work was done as one agent. Resources were only `PROBLEM.md`, `CONTRACT.md`,
`case.json`, local Python standard-library functionality, NumPy 2.4.6 and
SciPy 1.17.1. Tools used were workspace file reads/patches and local shell
execution. There was no network access, external source retrieval, other
model, package installation, or delegation. No particular published strip
crack solution was recognized or consulted; the named Dugdale and Maxwell
models were already supplied in the prompt. The analytical representation
was derived here from the stated mixed boundary problem.

The first recorded work timestamp was 2026-10-22 14:43:14 UTC. The final
solver run completed before 14:56:58 UTC; total work including packaging
was approximately 16 minutes (final delivery audit at 14:59:00 UTC), well
within the six-hour budget. The final
solver runtime was 41.5 seconds on the supplied environment. The final
generation was also run with Python warnings promoted to errors and
completed successfully; `verify.py` then passed all assertions.
