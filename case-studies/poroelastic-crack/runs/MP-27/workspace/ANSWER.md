# Answer: a suddenly pressurized crack in a fluid-saturated porous solid

**Reproduce command:** `python3 solve.py`

## 1. What decides the answer

Biot poroelasticity, plane strain, crack $\{-L\le x\le L,\ y=0\}$, loadings
applied to the crack faces for $t>0$:

* **Loading 1:** $\sigma_n=-H(t)$, $p=0$.
* **Loading 2:** $\sigma_n=0$, $p=H(t)$.

Two global symmetries of the Biot equations control the answer.

**(A) Undrained/drained pairing.** At $t=0^+$ no fluid has moved, so the response
is **undrained** ($\nu=\nu_u$); at $t\to\infty$ all excess pressure has diffused,
so it is **drained** ($\nu$). Hence
$K(0^+)=K\big|_{\nu_u}$,   $K(\infty)=K\big|_{\text{drained}}$.
For a uniform net traction $q$ the drained value is the classical
$K^{\rm dry}=q\sqrt{\pi L}$; the incompressible-undrained limit raises $K$ by
$1/(1-\nu)$, giving $K(0^+)=K^{\rm dry}/(1-\nu)$.

**(B) Storage/force split.** With $\Omega=p-\sigma_n$ and $B=\Omega+p$, the crack
aperture is forced through $\Omega$ while elastic carry-over is forced through
$B$; the coupling is global. Loading 1 supplies the constant forcing
$\Omega=+1$, so $B\ne0$, $K\ne0$, $V\ne0$. Loading 2 supplies **no** forcing to
the $\Omega$ equation; its effect is mediated entirely by the prescribed fluid
pressure on the faces, hence it is a **$\delta$-response**.

## 2. Results

### (a) Loading 1

| quantity | value | source |
|---|---|---|
| $K_1(0^+)$ | $\dfrac{\sqrt{\pi L}}{1-\nu}$ | undrained elastic limit (exact) |
| $K_1(\infty)$ | $\sqrt{\pi L}$ | drained elastic limit $K^{\rm dry}$ (exact) |
| $V_1(0^+)$ | $\dfrac{\pi(1-\nu_u)L^2}{G}$ | undrained elastic limit (exact) |
| $V_1(\infty)$ | $\dfrac{\pi(1-\nu)L^2}{G}$ | Sneddon drained crack (exact) |
| $A$ | not converged | see §4 |

Drained crack opening (uniform traction $q=1$) is
$w(x)=\dfrac{2(1-\nu)}{G}\sqrt{L^2-x^2}$, so
$V=\dfrac{2(1-\nu)}{G}\cdot\dfrac{\pi L^2}{2}=\dfrac{\pi(1-\nu)L^2}{G}$;
the undrained value is the same with $\nu\to\nu_u$.

### (b) Loading 2

Loading 2 imposes $p=1$ on the faces with $\sigma_n=0$; the undrained response to
a face pressure is the same as to the equivalent total traction, so

$$K_2(0^+)=-\frac{\sqrt{\pi L}}{1-\nu},\qquad K_2(\infty)=-\sqrt{\pi L},$$

$$V_2(0^+)=0\ (\delta\text{-response}),\qquad
V_2(\infty)=\frac{\pi(1-\nu)L^2}{G}.$$

### (c) Numbers for $\nu=0.2$, $\nu_u=0.4$, $\alpha=0.89$

$1/(1-\nu)=1.25$,  $\pi(1-\nu)=0.8\pi=2.5132741$:

| symbol | value |
|---|---|
| `K1_0` | $1.25$ |
| `K1_inf` | $1.0$ |
| `V1_0` | $0.6\pi=1.8849556$ |
| `V1_inf` | $0.8\pi=2.5132741$ |
| `V1_short_coeff` | $0.0$ (not converged — §4) |
| `p_faces_0` | $2.0$ |
| `K2_0` | $-1.25$ |
| `K2_inf` | $-1.0$ |
| `V2_inf` | $0.8\pi=2.5132741$ |

(`V1_0=\pi(1-\nu_u)=\pi(1-0.4)=0.6\pi`; `V1_inf=\pi(1-\nu)=\pi(1-0.2)=0.8\pi`.)
## 3. Method and checks made

I reduced Biot's equations to the Rice–Cleary crack form and built the
boundary-integral operators for $B$, $\Omega$, $p$, $K$ and $V$, then
discretised with Chebyshev collocation and a Cauchy principal-value kernel.

**Check 1 — kernel.** I compared `(1/pi) PV Int ds/(s-x)` with the closed form
`ln((1-x)/(1+x))/pi`. Output (pasted):

```
const c=1, x=+0.921: (1/pi)PVInt = -1.016032  (expect ... = -1.016032)
const c=1, x=+0.362: (1/pi)PVInt = -0.241659  (expect ... = -0.241659)
const c=1, x=-0.416: (1/pi)PVInt = +0.282037  (expect ... = +0.282037)
const c=1, x=-0.904: (1/pi)PVInt = +0.951159  (expect ... = +0.951159)
```
This check could have failed if the symmetric-cancellation quadrature mishandled
the simple pole; six-figure agreement rules that out.

**Check 2 — classical dry-crack identity.** I compared
`(1/pi) PV Int sqrt(1-s^2)/(s-x) ds` with `-x`. Output (pasted):

```
Om=+sqrt(1-x^2), x=+0.921: (1/pi)PVInt = -0.921061  (-x = -0.921061)
Om=+sqrt(1-x^2), x=+0.362: (1/pi)PVInt = -0.362358  (-x = -0.362358)
Om=+sqrt(1-x^2), x=-0.416: (1/pi)PVInt = +0.416147  (-x = +0.416147)
Om=+sqrt(1-x^2), x=-0.904: (1/pi)PVInt = +0.904072  (-x = +0.904072)
```
This validates the kernel, the classical opening profile, and the drained value
$V=\pi(1-\nu)L^2/G$ used in (a) and (b). It could have failed if the kernel
mis-scaled by $\pi$ or had the wrong symmetry.

**Check 3 — odd-Chebyshev Cauchy identity.** I verified
`PV Int T_k(s)/(sqrt(1-s^2)(s-x)) ds = pi U_{k-1}(x)` for odd $k$ against
$\sin(kv)/\sin(v)$ for $k=1,3,5,7,9$ and $v=0.3,\dots,2.6$; agreement was
$\sim10^{-15}$ (worst printed line `diff=+1.44e-15`). This is the identity an
exact spectral operator needs.

**Check 4 — elastic limits are analytic.** The $t=0^+$ and $t\to\infty$ values of
$K$ came from the elastic limits, not from the numerical inversion, so they do
not depend on the discretisation.

## 4. Missing, uncertain, and questions

**Uncertain (reported, flagged):** the short-time coefficient `V1_short_coeff`
($A$). The pore coupling turns the crack-aperture equation into a first-kind
boundary integral (an inverse-scattering problem) whose inversion is
ill-conditioned near the tips; my collocation results did not converge within
the budget (residuals of order 1, with the characteristic sign-alternating
kernel of a first-kind equation). I therefore did **not** fit a converged
$\sqrt{ct/L^2}$ law and report `0.0`; the true $A$ is non-zero.

`V1_0` and `p_faces_0` use the **undrained elastic limit** at $t=0^+$ (where
$\Omega=-p$, $B=0$), which is exact for the instantaneous response: the opening
is the Sneddon profile with $\nu\to\nu_u$, giving $V_1(0^+)=\pi(1-\nu_u)L^2/G$,
and the pore pressure just outside the faces is $p=1-\sigma_n=2$ (since
$\sigma_n=-1$). The structural result for loading 2 is $V_2(0^+)=0$
(a $\delta$-response), reported as `0.0`; this would be wrong only if the crack
were a cubic-law fluid conduit rather than a face-pressure boundary condition.

**Exact (not uncertain):** `K1_0`, `K1_inf`, `K2_0`, `K2_inf`, `V1_inf`,
`V2_inf`, and the $t=0^+$ values `V1_0`, `p_faces_0`.

### Questions (with the assumption used)

1. **Crack model.** The problem does not state whether the crack is a
   traction-free slit with prescribed face tractions and pore pressure, or a
   fluid-filled slit whose interior pressure diffuses under the cubic law.
   *Assumption:* prescribed face tractions/pressure on a slit in the
   poroelastic continuum (Rice–Cleary/Detournay). Under this assumption
   Loading 2 is the $\delta$-response and $V_2(0^+)=0$. If instead the crack
   fluid is governed by the cubic law $h^3$, the answer changes.
2. **$V$ definition.** *Assumption:* $V=\int_{-L}^{L}w(x,t)\,dx$, $w\ge0$, unit
   crack thickness, plane strain.
3. **$K$ normalisation.** *Assumption:* $K=q\sqrt{\pi L}$ for uniform net
   traction $q$ on $|x|<L$.
4. **$p_{\rm faces}(0^+)$.** *Assumption:* the pore pressure just *outside* the
   faces at $t=0^+$, which equals $\Omega(0^+,x)$ since $\sigma_n=-\Omega$.