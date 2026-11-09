# Reference solution

The answers are **$\Sigma=133.7302861\ \mathrm{MPa}$** and
**$L=4.085609138\ \mathrm{mm}$ at each end**. Grade each to ±0.3% using the
unrounded values in [reference_results.json](reference_results.json).

The calculation is implemented in [reference.py](reference.py). Run:

    python T1/reference.py

Only NumPy and SciPy are required; there are no data-file or network dependencies.

## Mechanical calculation

Use $E'=E/(1-\nu^2)=141.0590278\ \mathrm{GPa}$ and
$K_0=\sqrt{E'\Gamma_0}=1.761618180\ \mathrm{MPa}\sqrt{\mathrm m}$.
Let $c=a+L$ be the half-length including the connector zones and $\eta=c/a$.
The finite-crack stress intensity and opening at the original slit root are

\[
K_{\rm tip}=\sqrt{\pi c}\left[\Sigma-\frac{2\sigma_b}{\pi}
                  \cos^{-1}\!\left(\frac a c\right)\right],
\]

\[
\frac{E'\delta(a)}{\sigma_b a}
 =4\sqrt{\eta^2-1}\left[\frac{\Sigma}{\sigma_b}
                   -\frac2\pi\cos^{-1}(1/\eta)\right]
  +\frac8\pi\ln\eta .
\]

Quasistatic outer-front growth requires $K_{\rm tip}=K_0$, while first connector
failure requires $\delta(a)=\delta_c$. With
$u=\ln\eta$ and $k=K_0/(\sigma_b\sqrt{\pi a})$, eliminate $\Sigma$ to obtain

\[
4k\sqrt{2\sinh u}+\frac8\pi u=\frac{E'\delta_c}{\sigma_b a},
\qquad
\frac{\Sigma}{\sigma_b}=\frac2\pi\cos^{-1}(e^{-u})+ke^{-u/2}.
\]

The unique positive opening root is $u=0.994117006342$, giving
$c=6.485609138\ \mathrm{mm}$. The dimensionless load derivative
$d(\Sigma/\sigma_b)/du$ is $0.2150159321$ at this point. It is positive throughout
the preceding branch, so an earlier load maximum has not been mistaken for the
requested event.

The material length $E'\delta_c/\sigma_b$ is $7.934570313\ \mathrm{mm}$,
$a/(E'\delta_c/\sigma_b)=0.3024738462$, and $L/a=1.702337141$.
Replacing the coupled zone by the fully developed energy
$\Gamma_0+\sigma_b\delta_c=1462\ \mathrm{J/m^2}$ and a traction-free crack of
half-length $a$ would give $165.3840307\ \mathrm{MPa}$, an overprediction of
23.66984%. That approximation does not answer this finite-zone problem.

## Independent calculation and achieved agreement

The second route constructs a crack of half-length $c$ with face pressure

\[
q(s)=\begin{cases}
\Sigma,& |s|<a,\\
\Sigma-\sigma_b,& a<|s|<c.
\end{cases}
\]

It integrates the symmetric crack weight function for $K_{\rm tip}$ and the
exact crack-opening Green kernel for $\delta(a)$ by adaptive quadrature:

\[
K_{\rm tip}=2\sqrt{\frac c\pi}\int_0^c
                         \frac{q(s)}{\sqrt{c^2-s^2}}\,ds,
\]

\[
\delta(x)=\frac4{\pi E'}\int_{-c}^{c}q(s)
\ln\left|\frac{c^2-xs+\sqrt{(c^2-x^2)(c^2-s^2)}}{c(x-s)}\right|ds.
\]

A separate scalar root search imposes the opening cutoff, choosing $\Sigma$
from the numerically integrated tip condition. It does not call the closed
finite-crack opening formula above.

| Output | Closed-form reduction | Green/weight quadrature |
|---|---:|---:|
| $\Sigma$ (MPa) | 133.730286053811 | 133.730286053811 |
| $L$ (mm per end) | 4.085609138246 | 4.085609138246 |
| Root opening (µm) | 9.000000000000 | 9.000000000000 |

The observed relative differences are $2.2\times10^{-16}$ for stress and
$4.4\times10^{-16}$ for length. These are numerical agreement figures for the
specified ideal model, not an estimate of experimental accuracy.

Additional implemented checks:

- Uniform pressure gives the independently known elliptical opening
  $\delta(x)=4q\sqrt{c^2-x^2}/E'$; relative quadrature error is
  $1.6\times10^{-15}$. This checks the full-gap factor and kernel normalization.
- With $\Gamma_0=0$, the result reduces to
  $u=\pi E'\delta_c/(8\sigma_b a)$, agreeing to floating-point precision.
- For $\Gamma_0=0$ and $a=10^5 E'\delta_c/\sigma_b$, stress approaches
  $\sqrt{E'\sigma_b\delta_c/(\pi a)}$ within $6.6\times10^{-7}$ relative error.

All checks execute as assertions. Every reported reference quantity was checked.
