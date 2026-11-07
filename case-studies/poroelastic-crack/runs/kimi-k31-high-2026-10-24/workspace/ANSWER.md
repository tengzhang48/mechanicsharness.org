# Suddenly pressurized crack: derivation and results

Let the common loading magnitude be \(P\). The numerical request uses \(P=1\), and reports \(K/(P\sqrt{\pi L})\) and \(GV/(P L^2)\). Tension and opening are positive.

## 1. Elastic crack result

For a Griffith crack of length \(2L\) in a plane-strain elastic solid with Poisson ratio \(\nu_e\), a uniform opening traction \(P\) gives

\[
K=P\sqrt{\pi L},\qquad
\Delta u(x)=\frac{2(1-\nu_e)}{G}P\sqrt{L^2-x^2},
\]

where \(\Delta u=u_y(x,0^+)-u_y(x,0^-)\). Consequently,

\[
V=\int_{-L}^{L}\Delta u(x)\,dx
=\frac{\pi(1-\nu_e)}{G}P L^2.
\]

The instantaneous poroelastic response uses \(\nu_e=\nu_u\), whereas the fully drained response uses \(\nu_e=\nu\).

## 2. Loading 1: instantaneous response and face pore pressure

At \(t=0^+\), diffusion has not yet acted in the solid, so \(\zeta=0\). The fluid-content constitutive equation therefore gives

\[
p=-\beta\,\sigma_{kk},\qquad
\beta=\frac{\nu_u-\nu}{\alpha(1-2\nu)}.
\]

Substitution into the strain law eliminates \(p\):

\[
2G e_{ij}=\sigma_{ij}-\nu_u\sigma_{kk}\delta_{ij}.
\]

Thus the instantaneous solid is elastic with Poisson ratio \(\nu_u\). For the uniformly loaded Griffith crack, the stresses on the solid side of the crack faces are \(\sigma_{yy}=\sigma_{xx}=-P\), so \(\sigma_{kk}=-2P\). Hence

\[
\boxed{\frac{p_{\rm faces}(0^+)}{P}
=\frac{2(\nu_u-\nu)}{\alpha(1-2\nu)}}.
\]

This is the undrained pressure immediately inside the solid before drainage. The imposed value exactly on the hydraulically open face is \(p=0\) for \(t>0\); the ideal step loading creates the initial surface discontinuity that the diffusion layer subsequently resolves.

It follows that

\[
\boxed{K_1(0^+)=P\sqrt{\pi L}},\qquad
\boxed{V_1(0^+)=\frac{\pi(1-\nu_u)}{G}P L^2}.
\]

At long times the load-induced pore pressure has drained away, so the drained elastic constants apply:

\[
\boxed{K_1(\infty)=P\sqrt{\pi L}},\qquad
\boxed{V_1(\infty)=\frac{\pi(1-\nu)}{G}P L^2}.
\]

## 3. Loading 1: leading early-time growth

Away from the crack tips, the early diffusion field is one-dimensional in the coordinate \(y\ge0\) normal to either face. Write changes from the instantaneous state as \(\Delta(\cdot)\). The normal traction is fixed, so equilibrium gives \(\Delta\sigma_{yy}=0\). Uniformity along the face gives \(\Delta e_{xx}=0\). The \(xx\) strain equation therefore yields

\[
(1-\nu)\Delta\sigma_{xx}+\alpha(1-2\nu)\Delta p=0,
\]

and the \(yy\) strain equation then reduces to

\[
\Delta e_{yy}=D\,\Delta p,\qquad
D=\frac{\alpha(1-2\nu)}{2G(1-\nu)}.
\]

The stress term in the mass balance cannot be omitted. Using the above relation for \(\Delta\sigma_{xx}\), it becomes

\[
\frac{\alpha^2(1-2\nu)^2(1-\nu_u)}
{(\nu_u-\nu)(1-\nu)}\,\frac{\partial p}{\partial t}
=2G\kappa\,\frac{\partial^2p}{\partial y^2},
\]

which is exactly \(p_t=c p_{yy}\) for the diffusivity \(c\) given in the question.

Initially \(p=p_{\rm faces}(0^+)\equiv p_0\) in the solid, while \(p(0,t)=0\). Thus

\[
p(y,t)=p_0\,\operatorname{erf}\!\left(\frac{y}{2\sqrt{ct}}\right),
\qquad
\Delta p=-p_0\,\operatorname{erfc}\!\left(\frac{y}{2\sqrt{ct}}\right).
\]

Using \(\int_0^\infty\operatorname{erfc}(\eta)\,d\eta=1/\sqrt{\pi}\),

\[
\int_0^\infty\Delta p\,dy=-2p_0\sqrt{\frac{ct}{\pi}}.
\]

Since the far field is unchanged, the upper-face displacement change is

\[
\Delta u_y(0^+)=-\int_0^\infty\Delta e_{yy}\,dy
=2D p_0\sqrt{\frac{ct}{\pi}}.
\]

The lower face contributes equally, and integration over crack length \(2L\) gives

\[
V_1(t)-V_1(0^+)
=8L D p_0\sqrt{\frac{ct}{\pi}}+o(\sqrt{ct}).
\]

But

\[
D p_0=\frac{\nu_u-\nu}{G(1-\nu)}P,
\]

so

\[
\boxed{V_1(t)=\frac{\pi(1-\nu_u)}{G}P L^2
+\frac{8(\nu_u-\nu)}{G(1-\nu)}P L\sqrt{\frac{ct}{\pi}}
+o(\sqrt{ct})}.
\]

Equivalently, \(V_1\simeq V_1(0^+)+A\sqrt{ct/L^2}\) in normalized units, with

\[
\boxed{A=\frac{8(\nu_u-\nu)}{(1-\nu)\sqrt{\pi}}}.
\]

The neglected tip regions occupy only \(O(\sqrt{ct})\) of the crack length and do not change the displayed leading term.

## 4. Loading 2

At \(t=0^+\), the pore pressure has not entered the solid and the zero-stress, zero-displacement state satisfies all bulk equations and the total-traction condition. Therefore

\[
\boxed{K_2(0^+)=0}.
\]

In the drained long-time limit, the appropriate traction on the solid is expressed through the effective stress \(\sigma_{ij}+\alpha p\,\delta_{ij}\). On a crack face, \(\sigma_n=0\) and \(p=P\), so the effective normal stress is \(+\alpha P\). Since an opening load \(P_o\) corresponds to face normal stress \(\sigma_n=-P_o\), this is equivalent to \(P_o=-\alpha P\). The pore pressure in the surrounding bulk drains toward zero at fixed spatial position in the long-time diffusive limit. Hence

\[
\boxed{K_2(\infty)=-\alpha P\sqrt{\pi L}},
\qquad
\boxed{V_2(\infty)=-\frac{\alpha\pi(1-\nu)}{G}P L^2}.
\]

The negative signs denote crack closure relative to the unloaded configuration.

## 5. Numerical evaluation for \(\nu=0.2\), \(\nu_u=0.4\), \(\alpha=0.89\)

The following values are per unit loading magnitude. \(K\) is in units of \(\sqrt{\pi L}\), and \(V\) is in units of \(L^2/G\).

| Quantity | Value |
|---|---:|
| \(K_1(0^+)\) | 1 |
| \(K_1(\infty)\) | 1 |
| \(V_1(0^+)\) | 1.884955592154 |
| \(V_1(\infty)\) | 2.513274122872 |
| \(A\) | 1.128379167096 |
| \(p_{\rm faces}(0^+)\) | 0.749063670412 |
| \(K_2(0^+)\) | 0 |
| \(K_2(\infty)\) | -0.89 |
| \(V_2(\infty)\) | -2.236813969356 |

The retained script `verify.py` independently checks the undrained elimination, zero fluid content at the faces, the stress-coupled diffusion coefficient, the normal-strain relation, the erfc integral, the elliptical crack area, and all numerical relations. It also writes `results.json`. No unresolved uncertainties remain for the stated ideal boundary-value problem.
