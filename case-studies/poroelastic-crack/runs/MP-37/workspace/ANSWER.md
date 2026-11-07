# Suddenly pressurized poroelastic crack

## Answers

Write

\[
d=\nu_u-\nu,\quad a=\alpha(1-2\nu),\quad
b=\frac{a}{1-\nu},\quad \eta=\frac{d}{1-\nu_u}.
\]

The results below use the **true-tip** definition of the stress intensity:
first extract the stress singularity at each positive time, then let time tend
to zero. For this permeable crack, the opposite order of limits gives a
different instantaneous loading-1 intensity. Both are stated explicitly below.

The normalized answers for part (c) are:

| Quantity | Result |
| --- | ---: |
| `K1_0` (true tip) | 0.807593003 |
| `K1_inf` | 1 |
| `V1_0` | 1.884955592 |
| `V1_inf` | 2.513274123 |
| `V1_short_coeff` | 2.440219579 |
| `p_faces_0` (outer undrained value) | 0.7490636704 |
| `K2_0` | 0 |
| `K2_inf` | -0.33375 |
| `V2_inf` | -0.8388052385 |

Here K is in units of \(\sqrt{\pi L}\), and V in units of \(L^2/G\).
The extra digits for the numerically factored quantities are provided for
reproduction, not as a rigorous error certification. A conservative numerical
uncertainty is about \(10^{-6}\) in these two dimensionless quantities.

For loading 1,

\[
K_{1,0}^{\rm true}=\frac{n_{11}(\eta)}{1+\eta}\sqrt{\pi L},
\qquad K_{1,\infty}=\sqrt{\pi L},
\]
\[
V_{1,0}=\frac{\pi(1-\nu_u)L^2}{G},\qquad
V_{1,\infty}=\frac{\pi(1-\nu)L^2}{G},
\]
\[
V_1(t)=\frac{L^2}{G}\left[
\pi(1-\nu_u)+A\sqrt{\frac{ct}{L^2}}
 +o\left(\sqrt{\frac{ct}{L^2}}\right)\right],
\]
\[
\boxed{A=\frac{4d}{(1+\eta)\sqrt\pi}
\left[4+\pi\sqrt{\frac2\eta}\,n_{21}(\eta)\right]}.
\]

The functions \(n_{11},n_{21}\) are defined by an explicit canonical
matrix factorization below; no unestablished elementary approximation is used.
For the specified parameters,
\(n_{11}=1.0767906703323102\) and \(n_{21}=0.4169594526473146\).

The loading-1 bulk undrained face-adjacent pressure is

\[
p_{\rm faces,0}^{\rm outer}=\frac{2d}{a}.
\]

For loading 2,

\[
K_{2,0}=0,\qquad
K_{2,\infty}=-\frac{a}{2(1-\nu)}\sqrt{\pi L},\qquad
V_{2,\infty}=-\frac{\pi aL^2}{2G}.
\]

Negative K and V are signed responses of the prescribed linear traction model.
Loading 2 causes closure, not opening, when mechanical face traction is zero.

## Reproduction

Run from this folder:

```bash
OPENBLAS_NUM_THREADS=1 python3 solve.py
```

This recomputes the factorization and all reported numbers, writes
`results.json`, and prints its complete contents. It requires Python 3 and
NumPy. The independent finite-crack check additionally requires SciPy.
`OPENBLAS_NUM_THREADS=1 python3 solve.py --diagnostics` also prints the
factorization convergence and identity checks.

## Model And Method

I use the supplied two-index plane-strain equations literally, with zero shear
traction and symmetric mode I. Linear elasticity supplies the bulk endpoint
responses. A Fourier-Laplace half-space reduction supplies the coupled tip
boundary layer. Canonical Wiener-Hopf factorization determines its mechanical
and hydraulic amplitudes. Undrained elastic reciprocity then gives the complete
leading volume change, including both the smooth faces and the two tips.

This method is needed because a planar diffusion layer alone misses an
equally large tip contribution, and replacing the entire tip field by an
undrained elastic solution interchanges singular limits.

### Bulk Response

Let \(s=\sigma_{11}+\sigma_{22}\). Compatibility and equilibrium give

\[
\nabla^2 h=0,\quad h=s+bp,
\qquad
2G\zeta=ah+D_f p,\quad
D_f=\frac{a^2(1-\nu_u)}{d(1-\nu)}.
\]

Set \(\beta=a/D_f\), so \(b\beta=\eta\).
With \(\zeta=0\), \(p=-ds/a\), and the mechanical constitutive relation
reduces to elasticity with Poisson ratio \(\nu_u\).
For a uniformly pressured elastic slit,

\[
w(x)=\frac{2(1-\nu_*)}{G}\sqrt{L^2-x^2},\qquad
V=\frac{\pi(1-\nu_*)L^2}{G},\qquad K=\sqrt{\pi L}.
\]

Thus \(\nu_* = \nu_u\) initially in the outer bulk, and \(\nu_* = \nu\)
at long times for loading 1. The bulk elastic stress trace on either face is
\(s=-2\), giving the positive outer pressure \(2d/a\).

### Surface Symbol

Use Laplace variable \(z>0\), Fourier wavenumber \(\xi\),
\(k=|\xi|\), \(m=\sqrt{k^2+z/c}\), and opening \(w=2u_y(x,0^+)\).
Zero shear traction gives the exact surface relations

\[
\widehat S=-\frac{Gk}{2(1-\nu)}
\left[1+\eta\frac{m-k}{m+k}\right]\widehat w
-\frac{bk}{m+k}\widehat p_0,
\]
\[
\widehat Q=-m\widehat p_0+
\frac{\beta Gk(m-k)}{1-\nu}\widehat w,
\qquad Q=p_{,y}(x,0^+).
\]

On the crack, S and \(p_0\) are prescribed. On the intact symmetry line,
\(w=0\) and \(Q=0\). These latter conditions must not be replaced by
zero pressure outside the crack.

For example, an Airy representation is

\[
\Phi=A_0 e^{-ky}-\frac{(1+\eta)H}{2k}ye^{-ky}
-\frac{bP}{m^2-k^2}e^{-my},\quad
p=P e^{-my}-\beta H e^{-ky}.
\]

Here \(H\) is the surface amplitude of h, and \(P=p_0+\beta H\).
Enforcing \(\Phi_{,y}(0)=0\) and integrating the yy strain from infinity
gives \(H=-Gkw/(1-\nu)\), from which the displayed symbol follows.

### Canonical Tip Factor

Scale distances by \(\delta_z=\sqrt{c/z}\); in the next formulas
\(k=|\xi|\), \(m=\sqrt{1+\xi^2}\), and \(t_\xi=k/(m+k)\).
Define

\[
D=1+\eta-\frac{2\eta k^2}{m(m+k)},\quad
k_+=\sqrt{-i\xi},\ k_-=\sqrt{i\xi},\quad
m_+=\sqrt{1-i\xi},\ m_-=\sqrt{1+i\xi},
\]
\[
N(\xi)=
\begin{pmatrix}
D & \sqrt{2\eta}\,t_\xi/(m_+k_+)\\
\sqrt{2\eta}\,t_\xi/(m_-k_-) & 1
\end{pmatrix},\qquad N(0)=\operatorname{diag}(1+\eta,1).
\]

Square roots have positive real part. Factor

\[
N=N_+N_-,\qquad N_\pm(\infty)=I,
\]

with each factor and its inverse analytic in its corresponding upper/lower
half-plane. Define \(n_{11}=[N_-(0)]_{11}\),
\(n_{21}=[N_-(0)]_{21}\).
This specifies general expressions for all parameters in the problem; it is
a matrix factorization, not merely a factorization of \(\det N\).

For a tip at x=0, crack on x<0, use Fourier transform
\(\widehat f=\int f(x)e^{i\xi x}\,dx\).
The normalized singular cell has opening \(W\sim\sqrt{-X}\) far behind
the tip, with lower-supported opening and flux and upper-supported stress
and pressure. Its Wiener-Hopf pole residue is
\(\Gamma(3/2)N_-(0)(1,0)^T\).
At the inner tip,

\[
W\sim n_{11}\sqrt{-X},\qquad
P\sim B_p\sqrt R\cos(\theta/2),\quad
B_p=-\sqrt{2/\eta}\,n_{21}.
\]

The inner elastic field is drained while the outer opening is undrained,
giving \(K_{1,0}^{\rm true}/\sqrt{\pi L}=n_{11}/(1+\eta)\).
An independent Laplace-domain configurational integral, including the
pressure-gradient energy, gives

\[
n_{11}^2+n_{21}^2=1+\eta.
\]

Equivalently the outer elastic energy flux equals the inner drained elastic
flux plus \(\pi\kappa A_p^2/(4z)\), where
\(p\sim A_p\sqrt r\cos(\theta/2)\).
This identity both checks the factorization and explains why a purely elastic
modulus-ratio estimate does not suffice.

Numerically, `factor_check.py` maps the real frequency axis to a circle by
\(\xi=\rho\tan(\theta/2)\), solves
\(B=I-P_-[(N-I)B]\) by FFT projection, and normalizes
\(N_-(0)=B(\infty)B(0)^{-1}\).
`solve.py` uses \(\rho=0.01\), three fourfold grid refinements, and
extrapolation in \(N_{\rm grid}^{-1/2}\) to remove the frequency-zero cusp
error. The factor entries are not hard-coded.

### Short-Time Volume

Eliminating p from the supplied equations gives

\[
2Ge_{ij}=\sigma_{ij}-\nu_u s\delta_{ij}
 +\frac{2Gd}{a}\zeta\delta_{ij}.
\]

Reciprocity with the unit-pressure undrained elastic crack, whose trace is
\(s^a\), gives the exact identity

\[
V_1-V_{1,0}=\frac da\int_\Omega\zeta s^a\,dA.
\]

Away from the tips, the face layer has
\(p-p_{\rm outer}\simeq-p_{\rm faces,0}\operatorname{erfc}
(|y|/(2\sqrt{ct}))\).
Using both faces and \(\int_0^\infty\operatorname{erfc}(y/(2\sqrt{ct}))dy
=2\sqrt{ct}/\sqrt\pi\) gives

\[
A_{\rm faces}=\frac{16d}{(1+\eta)\sqrt\pi}.
\]

At each tip \(s^a\sim\sqrt{2L/r}\cos(\theta/2)\), so tip fluid-content
changes contribute at the same order. In the normalized tip cell put
\(H=-2kW\) on the surface and \(Z=P+H\). In the body,
\(\nabla^2H=0\), \((\nabla^2-1)Z=0\), \(\nabla^2P=Z\).
Fluid content is proportional to Z, so changes in h do not require a separate
extra bulk term. Green's identity with
\(f=R^{-1/2}\cos(\theta/2)\) yields

\[
I=\int_{\rm slit\ plane} Zf\,dA=-\pi B_p
=\pi\sqrt{2/\eta}\,n_{21}.
\]

The faces contribute zero to this Green integral, the outer-circle term
vanishes by undrained matching, and the inner-circle term is \(-\pi B_p\).
Including two tips and inverse-Laplace factor
\(\mathcal L^{-1}[z^{-3/2}]=2\sqrt t/\sqrt\pi\) gives
\(A_{\rm tips}=4d I/[(1+\eta)\sqrt\pi]\), and hence the reported A.

For small coupling,
\(n_{21}=\sqrt{\eta/2}+O(\eta^{3/2})\),
\(n_{11}=1+\eta/4+O(\eta^2)\). Therefore

\[
K_{1,0}^{\rm true}/\sqrt{\pi L}=1-3\eta/4+O(\eta^2),\qquad
A=\frac{16+4\pi}{\sqrt\pi}d+O(d^2)
\]

at fixed drained Poisson ratio. The weak-coupling approximation is not used
for part (c).

### Loading 2 Limits

Initially the imposed pressure occupies a layer of vanishing width. In its
tip cell the pressure stays bounded, although its gradient is singular.
Stress is O(1) on a length scale \(\sqrt{ct}\), so
\(K_2=O((ct)^{1/4})\to0\).

At long times there is no decaying harmonic pressure field with p=1 on a
finite slit in an infinite two-dimensional domain. Instead the pressure halo
expands; at each fixed position p tends to 1, nonuniformly in space.
A pressure particular solution can be written

\[
u_i^p=\frac b{2G}\psi_{,i},\quad \nabla^2\psi=p,
\qquad \sigma_{ij}^p=b(\psi_{,ij}-p\delta_{ij}).
\]

The asymptotically radial halo selects
\(\psi_{,ij}\to\delta_{ij}/2\), hence a smooth local total stress
\(\sigma_{ij}^p\to-b\delta_{ij}/2\).
The correction enforcing zero face traction is the drained elastic slit
under signed opening pressure \(-b/2\). Its displacement jump supplies
\(V_{2,\infty}\); the smooth particular displacement has no jump.

Taking \(m\to k\) in the nonzero-frequency surface symbol gives the same
traction term \(-bp_0/2\), but applying that limit directly to a constant
zero-frequency pressure would be invalid. Matching through the expanding
logarithmic halo supplies the stated \(-b/2\) background.

## Questions

- Shear traction was not explicitly prescribed. I assumed zero on both faces,
  with displacement and pressure symmetry across y=0 outside the slit.
- `K1_0` is ambiguous without an order of limits. I used
  \(\lim_{t\downarrow0}\lim_{r\downarrow0}\sqrt{2\pi r}
  \sigma_{yy}(L+r,0,t)\). If instead the instantaneous outer undrained
  elastic field is intended, use **`K1_0 = 1`**. That is
  \(\lim_{r\downarrow0}\lim_{t\downarrow0}\), not the true-tip value
  at positive time. The volume endpoints are unaffected by this distinction.
- The actual loading-1 face pressure trace is exactly zero. The requested
  `p_faces_0` was interpreted as the initially induced solid pressure before
  diffusion, outside the vanishing layer. If a literal boundary trace is
  intended, its value is **0** instead of the reported outer value.
- Loading 2 has no classical decaying two-dimensional steady-pressure solution.
  I used its transient long-time limit, retaining the far-field constraint.
- No contact/nonpenetration law is supplied. I allowed signed opening and
  negative volume as required by the stated linear model; enforcing contact
  would change loading 2.
- The permitted endpoint \(\nu_u=1/2\) is obtained by continuity. The
  numerical factorization used for (c) is away from that endpoint. General
  formulas are factorization-defined, rather than elementary closed forms.

## Checks And Uncertainty

Every numerical check is documented with the compared quantities, exact
pasted terminal output, and ways it could fail. See `CHECKS.md` for the primary
factorization and reproduction checks, and `finite_check.md` for the full
independent finite-crack audit, including failed and timed-out refinements.

The final extrapolated factor norm error is \(-4.67\times10^{-10}\).
Linear versus quadratic extrapolations differ by about \(1.4\times10^{-9}\)
in K and \(2.3\times10^{-9}\) in A, but this agreement alone is not a
rigorous discretization bound. Varying the frequency-map scale checks a
different source of numerical error.

The independent finite-crack mixed-boundary solve gives limiting estimates
around 2.438 to 2.441 before box extrapolation, consistent with the reported
2.440220 and inconsistent with the face-only 1.354055. Its finite-grid,
finite-Laplace-parameter, periodic-box, and empirical extrapolation errors are
larger than the primary factorization error. It validates the size and inclusion
of the tip term, not all displayed digits. No direct time-stepping or numerical
loading-2 long-time simulation was performed. No rigorous numerical error bound
or elementary factorization has been established.
