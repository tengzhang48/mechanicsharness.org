# Solution

Take each imposed loading to have unit magnitude. A different loading magnitude multiplies all reported stress-intensity factors, volumes, and induced pressures linearly. Let

\[
S=\sigma_{xx}+\sigma_{yy},\qquad a=\alpha(1-2\nu),\qquad d=\nu_u-\nu>0,
\]

so that the fluid-content equation is \(2G\zeta=aS+(a^2/d)p\).

## 1. Elastic crack facts and the two instantaneous limits

For an ordinary plane-strain elastic Griffith crack under a uniform compressive face traction \(\sigma_n=-q\), with elastic Poisson ratio \(\nu_e\),

\[
K^{(e)}=q\sqrt{\pi L},
\]

\[
h^{(e)}(x)=\frac{2(1-\nu_e)q}{G}\sqrt{L^2-x^2},
\qquad
V^{(e)}=\int_{-L}^{L}h^{(e)}(x)\,dx
=\frac{\pi(1-\nu_e)qL^2}{G}.
\]

Here \(h=u_y^+-u_y^-\) is the full opening. The stress field for this traction problem also gives

\[
S=-2q\quad\text{on the crack faces}.
\]

At \(t=0^+\), diffusion has not yet changed the fluid content in the solid, so \(\zeta=0\). Consequently,

\[
p^{(0)}=-\frac{d}{a}S,
\]

and the constitutive law becomes

\[
2G e_{ij}=\sigma_{ij}-\nu_uS\delta_{ij}.
\]

Thus the instantaneous mechanical response is elastic with \(\nu_e=\nu_u\). At long time under loading 1, the pore-pressure disturbance has drained away and the response is elastic with \(\nu_e=\nu\).

### Loading 1: endpoints and initial face pressure

Using \(q=1\) in the elastic formulas gives

\[
\boxed{K_1(0^+)=K_1(\infty)=\sqrt{\pi L}},
\]

\[
\boxed{V_1(0^+)=\frac{\pi(1-\nu_u)L^2}{G}},\qquad
\boxed{V_1(\infty)=\frac{\pi(1-\nu)L^2}{G}}.
\]

Immediately next to either face, before diffusion, \(S=-2\). Therefore the induced pore pressure is positive and equals

\[
\boxed{p_{\rm faces}(0^+)=\frac{2(\nu_u-\nu)}{\alpha(1-2\nu)}}.
\]

This is the value just inside the solid; the imposed boundary value at the drained face itself is zero for \(t>0\).

## 2. Leading short-time growth of \(V_1\)

Compatibility, equilibrium, and the constitutive law give

\[
(1-\nu)\nabla^2S+a\nabla^2p=0.
\]

Hence \(S+\eta p\), with \(\eta=a/(1-\nu)\), is harmonic. In the thin diffusion layer next to a flat portion of either crack face, let \(y\geq0\) be the inward normal coordinate. To leading order, the harmonic combination is constant across the layer:

\[
S+\eta p=H.
\]

Define

\[
b=\frac{a^2}{d},\qquad
B=b-a\eta=\frac{a^2(1-\nu_u)}{d(1-\nu)}.
\]

Then \(2G\zeta=aS+b p=aH+Bp\). Mass balance therefore reduces locally to

\[
B p_t=2G\kappa p_{yy},
\]

or \(p_t=c p_{yy}\), because

\[
\frac{2G\kappa}{B}
=\frac{2\kappa G(1-\nu)(\nu_u-\nu)}
{\alpha^2(1-2\nu)^2(1-\nu_u)}=c.
\]

The initial value just inside the face is \(p_f=2d/a\), while the boundary value is zero. The leading boundary-layer solution is consequently

\[
p(y,t)=p_f\,\operatorname{erf}\!\left(\frac{y}{2\sqrt{ct}}\right),
\qquad
p-p_f=-p_f\,\operatorname{erfc}\!\left(\frac{y}{2\sqrt{ct}}\right),
\]

and

\[
\zeta=\frac{B}{2G}(p-p_f)
=-\frac{Bp_f}{2G}\operatorname{erfc}\!\left(\frac{y}{2\sqrt{ct}}\right).
\]

To convert this fluid-content change into crack volume, use an undrained elastic reciprocal identity. Define

\[
\Lambda=\frac{2Gd}{a(1-2\nu_u)},\qquad
\tau_{ij}=\sigma_{ij}+\Lambda\zeta\delta_{ij}.
\]

The constitutive law can then be written as the undrained elastic law

\[
2G e_{ij}=\tau_{ij}-\nu_u\tau_{kk}\delta_{ij}.
\]

Betti's theorem, using the unit-pressure undrained elastic crack solution \(u^{(u)}\) as the auxiliary solution, gives

\[
V_1=V_1(0^+)+\Lambda\int_\Omega \zeta\,\theta_u\,dA,
\qquad \theta_u=\nabla\cdot u^{(u)}.
\]

On the crack faces,

\[
\theta_u=\frac{(1-2\nu_u)(-2)}{2G}
=-\frac{1-2\nu_u}{G},
\qquad
\Lambda\theta_u=-\frac{2d}{a}.
\]

Thus, in either boundary layer,

\[
\Lambda\theta_u\zeta
=\frac{2d(1-\nu_u)}{G(1-\nu)}
\operatorname{erfc}\!\left(\frac{y}{2\sqrt{ct}}\right).
\]

There are two faces, the loaded crack length is \(2L\), and

\[
\int_0^\infty \operatorname{erfc}\!\left(\frac{y}{2\sqrt{ct}}\right)dy
=\frac{2\sqrt{ct}}{\sqrt\pi}.
\]

Integration therefore gives the leading term

\[
V_1(t)-V_1(0^+)
\simeq \frac{16(\nu_u-\nu)(1-\nu_u)}{G(1-\nu)}
L\sqrt{\frac{ct}{\pi}}.
\]

Equivalently, in the normalization requested in the problem,

\[
\boxed{
\frac{GV_1}{L^2}
=\pi(1-\nu_u)+A\sqrt{\frac{ct}{L^2}}
+o\!\left(\sqrt{\frac{ct}{L^2}}\right)
}
\]

with

\[
\boxed{
A=\frac{16(\nu_u-\nu)(1-\nu_u)}{(1-\nu)\sqrt\pi}
}.
\]

The reciprocal calculation includes both the displacement accumulated across the diffusion layers and the accompanying undrained outer response.

## 3. Loading 2

Initially there is no total traction and no bulk deformation, so the zero-thickness pressure boundary layer has no integrated mechanical effect:

\[
\boxed{K_2(0^+)=0}.
\]

For the long-time drained response, use the effective stress

\[
\bar\sigma_{ij}=\sigma_{ij}+\alpha p\delta_{ij}.
\]

The sum of loadings 1 and 2 has face conditions \(\sigma_n=-1\), \(p=1\), hence effective normal traction

\[
\bar\sigma_n=-1+\alpha=-(1-\alpha).
\]

It is therefore a drained elastic crack under effective pressure magnitude \(1-\alpha\), giving

\[
K_{1+2}(\infty)=(1-\alpha)\sqrt{\pi L},
\qquad
V_{1+2}(\infty)=(1-\alpha)\frac{\pi(1-\nu)L^2}{G}.
\]

Subtracting the loading-1 results yields

\[
\boxed{K_2(\infty)=-\alpha\sqrt{\pi L}},
\]

\[
\boxed{V_2(\infty)=-\alpha\frac{\pi(1-\nu)L^2}{G}}.
\]

The negative signs are the formal signs for this auxiliary loading: positive face pore pressure with zero total normal traction produces a negative mode-I factor and negative opening volume.

## 4. Numerical evaluation for \(\nu=0.2\), \(\nu_u=0.4\), \(\alpha=0.89\)

The numerical fields below use the requested normalizations: \(K\) is divided by \(\sqrt{\pi L}\), and \(V\) is divided by \(L^2/G\). The induced pressure is per unit loading magnitude.

| Quantity | Value |
|---|---:|
| `K1_0` | 1.0 |
| `K1_inf` | 1.0 |
| `V1_0` | 1.884955592153876 |
| `V1_inf` | 2.5132741228718345 |
| `V1_short_coeff` | 1.354055000514615 |
| `p_faces_0` | 0.7490636704119851 |
| `K2_0` | 0.0 |
| `K2_inf` | -0.89 |
| `V2_inf` | -2.236813969355933 |

In particular, the short-time expression is

\[
\frac{GV_1}{L^2}\simeq1.884955592153876
+1.354055000514615\sqrt{\frac{ct}{L^2}}.
\]

## 5. Checks

- Setting \(\zeta=0\) in the original constitutive equations gives exactly the elastic law with \(\nu_u\); setting \(p=0\) gives the drained elastic law with \(\nu\).
- Integrating the elliptical elastic opening gives the two volume limits, and the near-tip stress field gives \(K=q\sqrt{\pi L}\).
- The boundary-layer pressure profile satisfies \(p_t=c p_{yy}\), with \(c\) identical to the diffusivity supplied in the statement.
- The reciprocal identity is consistent with the endpoint volume change. For the unit-pressure Griffith-crack stress field, \(\int_\Omega S^2\,dA=2\pi L^2\). At long time under loading 1, \(\zeta=aS/(2G)\), so the identity gives \(V_1(\infty)-V_1(0^+)=\pi dL^2/G\), as required.
- Loading superposition gives the numerical checks \(K_{1,\infty}+K_{2,\infty}=1-\alpha=0.11\) and \(V_{1,\infty}+V_{2,\infty}=(1-\alpha)\pi(1-\nu)=0.2764601535159015\) in the reported units.
- `verify.py` retains the constitutive, diffusivity, elastic-crack, reciprocal-identity, quadrature, short-time, superposition, and JSON-value checks. Running `python verify.py` passes all checks.

There are no unresolved uncertainties. The short-time result is a leading asymptotic result only; no claim is made here about the intermediate-time response.
