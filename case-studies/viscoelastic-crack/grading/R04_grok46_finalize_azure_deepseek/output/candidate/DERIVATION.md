# Derivation notes (working)

## Constitutive form
Standard linear solid: equilibrium spring $G_\infty$ parallel to a Maxwell branch
(spring $G_1=G_0-G_\infty$, dashpot $\eta$, $t_M=\eta/G_1$).
$$(1+t_M\partial_t)\tau_i=(G_\infty+t_MG_0\partial_t)\gamma_i .$$
Creep compliance $C(t)=\frac1{G_\infty}[1-(1-\lambda)e^{-t/t_c}]$, $\lambda=G_\infty/G_0$,
$t_c=t_M/\lambda$; free energy as given.

## Steady moving frame
$x=x_{lab}-vt$, stationary fields, material derivative $\dot{}=-v\partial_x$.
Eliminating the overstress gives
$$(G_\infty - v t_M G_0 \partial_x)\nabla^2 w = 0 .$$
With decay ahead ($\nabla^2 w\to0$ as $x\to+\infty$) this forces
$$\boxed{\nabla^2 w=0}$$ everywhere. So $w$ is harmonic; viscoelasticity enters only
through the traction/opening relation on $y=0$.

## Crack-plane operator
Half strip $0<y<h$, $\hat w(q,y)=\hat u(q)\sinh(q(h-y))/\sinh(qh)$, so
$\hat w_{,y}(0)=-q\coth(qh)\hat u$.
With $U=G_\infty w(x,0^+)/(\tau_c h)$, $X=\pi x/h$, $K_0=1/(\lambda V)$, $V=\pi vt_c/h$:
$$\hat T_2(Q)=a\,2\pi\delta(Q)-\ell(Q)\hat U(Q),\qquad
\ell(Q)=\frac{\pi}{\lambda}\,Q\coth(\pi Q)\,\frac{iQ-K_0\lambda}{iQ-K_0}.$$
Limits: $K_0\to\infty\Rightarrow \ell=\pi Q\coth(\pi Q)$ (relaxed, $G_\infty$);
$K_0\to0\Rightarrow \ell=\frac\pi\lambda Q\coth(\pi Q)$ (glassy, $G_0$).

## Local (bulk) form used numerically
$$\nabla^2 W=0,\qquad W(X,\pi)=a,\quad X=\pi x/h,\ Y=\pi y/h,$$
$$T_2=\pi W_{,Y}+S_2,\quad T_1=\pi W_{,X}+S_1,$$
$$(\partial_X-K_0)S_1=\tfrac{\pi(1-\lambda)}{\lambda}W_{,XX},\qquad
(\partial_X-K_0)S_2=\tfrac{\pi(1-\lambda)}{\lambda}W_{,XY},$$
with $S_i\to0$ as $X\to+\infty$, i.e.
$S_i(X,Y)=-\tfrac{\pi(1-\lambda)}{\lambda}\int_0^\infty W_{,Xi}(X+u,Y)e^{-K_0u}\,du$.
Boundary data on $Y=0$: $W=0$ ($X\ge0$), $T_2=1$ ($-\epsilon<X<0$), $T_2=0$ ($X<-\epsilon$);
far field $W\to(a/\pi)Y$ (ahead), $W\to a$ (behind).
Unknowns $a,\epsilon=\pi l_c/h$ from: bounded $T_2$ at $X=0$ and $2U(-\epsilon)=\bar\delta_c$.