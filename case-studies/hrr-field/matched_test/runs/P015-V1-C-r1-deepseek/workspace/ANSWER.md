# Answer: Slip and strain contours near a stationary slit

## Requested numbers (mm, MPa)

1. **Relative face slip** \(\Delta w(r)=w(r,+\pi)-w(r,-\pi)\):
   - \(r=0.0200\) mm: \(\Delta w = 1.3181\times10^{-3}\) mm
   - \(r=0.0800\) mm: \(\Delta w = 1.7937\times10^{-3}\) mm

2. **Signed stresses** at \(r=0.0200\) mm, \(\theta=60^\circ\):
   - \(\tau_{rz} = +38.59\) MPa
   - \(\tau_{\theta z} = +138.01\) MPa

3. **Radii where \(|\gamma|=0.0100\)**:
   - \(\theta=0^\circ\): \(r = 0.046894\) mm
   - \(\theta=90^\circ\): \(r = 0.025066\) mm
   - \(\theta=180^\circ\) (upper face, limit from material): \(r = 0.013398\) mm

(Computed values: \(A = 4.16524\times10^{-3}\); \(\Delta w(0.02)=1.318126778\times10^{-3}\);
\(\Delta w(0.08)=1.793693738\times10^{-3}\); \(\tau_{rz}=38.5866\); \(\tau_{\theta z}=138.0106\);
\(r_0=0.0468944\); \(r_{90}=0.0250661\); \(r_{180}=0.0133984\).)

---

## Governing model and why

The loading is antiplane shear (mode III): the only displacement is \(w(x,y)=u_z\), the slit lies
on \(\theta=\pm\pi\), and the free faces and odd loading fix a Mode-III symmetry. With
\(\boldsymbol\gamma=(w_{,x},w_{,y})=\nabla w\) and the given law
\(\boldsymbol\gamma=\gamma_0(\tau/\tau_0)^{n-1}\boldsymbol\tau/\tau_0\), inverting gives
\(\boldsymbol\tau = \tau_0\,(|\nabla w|/\gamma_0)^{1/n}\,\nabla w/|\nabla w|\). Equilibrium
\(\nabla\!\cdot\!\boldsymbol\tau=0\) becomes

\[
\nabla\cdot\bigl(|\nabla w|^{\,1/n-1}\nabla w\bigr)=0 ,
\]

i.e. \(w\) is \(p\)-harmonic with \(p=1+1/n=(n+1)/n=9/7\).

The near-tip leading field is the separable solution \(w=A\,r^{a}F(\theta)\). Homogeneity of the
\(p\)-Laplacian forces the **eigenvalue** \(a=(p-1)/p=1/(n+1)=2/9\), so that
\(\tau\sim r^{-1/(n+1)}\), \(\gamma\sim r^{-n/(n+1)}\) and the \(J\)-integrand is \(r\)-independent.
Using \(u=F'/(aF)\), the angular problem has the closed form (with \(a=1/(n+1)\))

\[
\theta(u)=\pi-\arctan u-\arctan(u/n),\qquad
F(u)=(1+u^{2})^{-a/2}\,(u^{2}+n^{2})^{-an/2},\qquad F'=a u F,
\]

which satisfies the traction-free faces \(F'(\pm\pi)=0\), the odd symmetry \(F(-\theta)=-F(\theta)\),
and the eigenvalue condition \(\int_0^\infty |d\theta/du|\,du=\pi\). (For \(n=1\) it reduces to the
classic linear Mode-III field \(F\propto\sin(\theta/2)\).)

The stored work for a proportional path is \(W=\int\tau\,d\gamma=\tfrac{n}{n+1}\tau\gamma\).

### How the loading fixes the magnitude

The amplitude \(A\) is the only free constant. The \(J\)-integral computed on a circle of any
radius \(r\) is \(r\)-independent and equals

\[
J=\tau_0\,\gamma_0^{-1/n}\,A^{(n+1)/n}\,D ,
\]

where \(D\) is a pure number from the angular field. Setting \(J=0.120\) N/mm gives
\(A=(0.120/J_{A=1})^{n/(n+1)}=4.16524\times10^{-3}\) mm\(^{\,n/(n+1)}\) (numerically
\(J_{A=1}=137.925\) N/mm).

---

## Assumptions changed and why

- **No elastic term**: used the specified pure power law exactly (\(W=\tfrac{n}{n+1}\tau\gamma\)).
- **Leading (asymptotic) field only**: geometry outside \(r=0.30\) mm was declared negligible, so
  the homogeneous separable near-tip field is used with no higher-order corrections.
- **Parity**: \(w\) is *odd* in \(\theta\) (from \(w(x,-y)=-w(x,y)\)), consistent with
  \(w\sim r^{a}\sin(\cdot)\); the stress function \(\phi\) is even. The sign of \(\tau_{\theta z}\)
  follows the stated positive basis \(e_\theta=(-\sin\theta,\cos\theta)\).

Nothing else was altered.

---

## Verification

- **Energy-flux / path independence**: the same \(J\)-integral integrated on circles
  \(r=0.005,0.02,0.1,0.3\) mm gives \(J=0.120000\) N/mm for all radii (with the fitted \(A\)).
- **Equilibrium residual**: \(\nabla\!\cdot\!\boldsymbol\tau\) computed by central differences is
  \(\lesssim 10^{-6}\) MPa/mm against \(|\boldsymbol\tau|\sim 100\) MPa (relative \(\sim10^{-9}\)).
- **Boundary condition**: \(\tau_{\theta z}\to 0\) on \(\theta=\pm\pi\) (face traction free).
- **Independent ODE integration**: integrating the angular \(p\)-Laplace ODE
  \(F''=-F\,\dfrac{c\,a^2F^2+(c+(p-2)a^2)F'^2}{a^2F^2+(p-1)F'^2}\), \(c=n/(n+1)^2\), with
  \(F(0)=0,F'(0)=1\) gives a shape proportional to the closed form (ratio \(1/a=4.5\)) and
  \(F'(\pi)\approx-3\times10^{-17}\), confirming the closed form and eigenvalue.
- **Work law**: \(\int_0^\gamma\tau\,d\gamma\) matches \(\tfrac{n}{n+1}\tau\gamma\) to \(10^{-12}\).

---

## Missing / uncertain

- The model returns the *local* field for a *given* \(J\); it does not know the far-field load, and
  the small-strain, stationary, monotonic assumptions are taken as given.
- **Do these results alone determine the load at which the slit starts to advance? No.** The
  results give the near-tip stress/strain/slip as functions of \(J=0.120\) N/mm, but the load
  corresponding to that \(J\) requires the global boundary-value solution, and the onset of advance
  additionally requires a fracture criterion (e.g. a critical \(J\) or critical crack-tip value),
  neither of which is supplied — so the initiation load is not determined by these local fields
  alone.
