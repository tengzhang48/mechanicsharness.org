# Steady crack growth in a viscoelastic strip

Determine the steady propagation states of the following crack problem. Derive the solution from the stated model and supply executable numerical results. The analytical or numerical method is yours to choose.

## Geometry and loading

A homogeneous strip occupies \(-\infty<x_{\rm lab}<\infty\), \(-h<y<h\), and is invariant in the third direction. Deformation is small and antiplane: the only displacement is \(w\) in the third direction. Inertia and body forces are absent. Let

\[
\gamma_1=w_{,x},\quad \gamma_2=w_{,y},\qquad
\tau_1=\tau_{x3},\quad\tau_2=\tau_{y3},\qquad
\tau_{1,x}+\tau_{2,y}=0.
\]

The crack propagates to the right at constant speed \(v>0\). Use the moving coordinate \(x=x_{\rm lab}-vt\), choosing its origin at the front of the cohesive zone. The intact interface is \(y=0,\ x\geq0\); the cohesive zone is \(-l_c<x<0\); the fully separated crack is \(x<-l_c\). The length \(l_c\) is unknown. Fields are stationary in this moving frame.

The upper and lower edges have fixed displacements \(w(x,\pm h)=\pm\Delta\), with \(\Delta>0\) unknown. Displacement is odd in \(y\). In the upper half strip,

\[
w(x\geq0,0)=0,\qquad
\tau_2(-l_c<x<0,0^+)=\tau_c,\qquad
\tau_2(x<-l_c,0^+)=0.
\]

Require bounded stress at the front of the cohesive zone. Material far ahead is fully relaxed: \(w\to(\Delta/h)y\), \(\gamma_1\to0\), \(\gamma_2\to\Delta/h\), \(\tau_1\to0\), and \(\tau_2\to G_\infty\Delta/h\) as \(x\to+\infty\). Far behind, the upper half has \(w\to\Delta\), zero stress and zero strain.

## Bulk material and fracture law

Each shear component obeys linear Boltzmann superposition with the shear creep compliance

\[
C(t)=\frac{1}{G_\infty}
\left[1-(1-\lambda)e^{-t/t_c}\right],\qquad t\geq0,
\qquad\lambda=G_\infty/G_0.
\]

Here \(G_0\) and \(G_\infty\) are the instantaneous and relaxed shear moduli, and \(t_c\) is the creep time appearing in this formula. Histories begin in the fully relaxed far-ahead state specified above. The material has an equilibrium spring \(G_\infty\) in parallel with a Maxwell branch whose spring modulus is \(G_0-G_\infty\). For calculating local energy loss, use its free energy density

\[
\psi=\frac{G_\infty}{2}(\gamma_1^2+\gamma_2^2)
+\frac{(\tau_1-G_\infty\gamma_1)^2+
(\tau_2-G_\infty\gamma_2)^2}{2(G_0-G_\infty)}.
\]

The interface obeys a rate-independent Dugdale law. There is no displacement jump before the positive shear traction reaches \(\tau_c\). During separation the traction remains \(\tau_c\), until the relative displacement reaches \(\delta_c\), when the traction drops to zero. Thus the steady solution must satisfy \(2w(-l_c,0^+)=\delta_c\). Intrinsic fracture energy is \(\Gamma_0=\tau_c\delta_c\).

## Requested solution

Use \(\lambda=0.1\) and \(\bar\delta_c=G_\infty\delta_c/(\tau_ch)=0.2\). Define

\[
X=\pi x/h,\quad Y=\pi y/h,\quad V=\pi vt_c/h,\quad
\epsilon=\pi l_c/h,\quad a=G_\infty\Delta/(\tau_ch).
\]

1. Derive the relation between loading and speed, together with the cohesive-zone length and stress and strain fields. Determine the range of positive speeds admitting steady solutions; report any finite endpoint and the limiting behavior there. State the slow-speed limit.
2. Compute \(\epsilon\), \(a\), and \(R=\Gamma/\Gamma_0\) at \(V=0.001,0.1,10\), where \(\Gamma\) is the loss of stored energy per unit advance between the relaxed far-ahead and stress-free far-behind states. Report normalized field samples at the coordinates in `case.json`, with \(T_i=\tau_i/\tau_c\) and \(E_i=G_\infty\gamma_i/\tau_c\).
3. Compute the bulk energy dissipated per unit crack advance by spatial integration of \(\phi_D=\tau_i\dot\gamma_i-\dot\psi\), where dots follow material points and repeated \(i\) is summed. Report \(\Gamma_D/\Gamma_0\), the integration bounds, refinement and tail checks, and the energy-balance discrepancy. An energy-loss value obtained only by subtracting \(\Gamma_0\) from \(\Gamma\) does not meet this integration requirement.
4. Supply numerical checks of equilibrium, the material law, interface and outer boundary conditions, and the critical separation. Identify which errors each check could detect. Report unresolved quantities explicitly.

Write `ANSWER.md`, runnable Python code and `results.json` as specified in `CONTRACT.md`. All physical inputs are above. `case.json` specifies numerical sampling only. No particular discretization, mapping or solution representation is required.
