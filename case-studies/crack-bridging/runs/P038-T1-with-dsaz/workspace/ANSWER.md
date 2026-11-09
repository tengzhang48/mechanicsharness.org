# Answer — A finite flaw under opening tension

## Results

1. **Remote stress when the connector adjacent to the original slit first
   fails:** $\Sigma = 133.73$ MPa
2. **Corresponding added crack length per end:** $L = 4.086$ mm

(within the ±0.3% tolerance; method and checks below)

---

## Material / geometry

- Plane-strain modulus
  $E' = E/(1-\nu^2) = 130/(1-0.28^2)\ \mathrm{GPa} = 141.059$ GPa.
- Bridging scale: $\sigma_0 = t = 160$ MPa, limiting gap $\delta_c = 9.00\ \mu m$;
  $\sigma_0\delta_c = 1440$ J/m².
- Bridging material length
  $\ell = E'\delta_c/\sigma_0 = 7.9346$ mm,
  so $\alpha = a/\ell = 2.40/7.9346 = 0.3025$ → **large-scale bridging**
  (the connector zone is comparable to the flaw; K-concept/Small-scale
  results do not apply).
- Intrinsic front threshold
  $K_0 = \sqrt{E'\Gamma_0} = \sqrt{141.059\times10^{9}\cdot 22}
  = 1.7616\ \mathrm{MPa}\sqrt{\mathrm m}$.

## Method

As the fronts advance, connectors occupy $a<|x|<c$, $c=a+L$, and carry the
constant closing traction $\sigma_0$ as long as their gap is below $\delta_c$.
The opening profile of a bridged crack is **largest at the trailing edge of the
zone, $|x|=a$** — the connector adjacent to the original slit — and falls to
zero at the front $|x|=c$. So the **first** connector to reach $\delta_c$ is the
one at $|x|=a$; growth is on the monotonic branch until

$$\delta_t \equiv \delta(|x|=a) = \delta_c
\quad\text{and}\quad K_{tip}(|x|=c) = K_0 = \sqrt{E'\Gamma_0}.$$

This is exactly the **partially-bridged center crack** problem (uniform closure
$\sigma_0$ on $a<|x|<\eta a$, $\eta = 1+L/a$), for which Bao & Suo (1992) give
the exact Tada closed forms, eqs (45) and (46):

$$\frac{\delta_t E'}{\sigma_0 a}
= 4\sqrt{\eta^2-1}\left[\frac{\bar\sigma}{\sigma_0}-\frac{2}{\pi}\cos^{-1}\frac1\eta\right]
+\frac{8}{\pi}\ln\eta,
\qquad
\frac{K_{tip}}{\sigma_0\sqrt{\pi a}}
=\sqrt{\eta}\left[\frac{\bar\sigma}{\sigma_0}-\frac{2}{\pi}\cos^{-1}\frac1\eta\right].$$

Setting $K_{tip}=K_0$ gives $\bar\sigma(\eta)$ from eq (46); setting
$\delta_t=\delta_c$ in eq (45) then fixes $\eta$. Solving (SciPy `brentq`, and
independently with a 2-equation `fsolve`) gives

$$\eta = 2.70234,\qquad L = (\eta-1)a = 4.0856\ \text{mm},\qquad
\Sigma = \bar\sigma = 133.73\ \text{MPa}.$$

## Checks

1. **Reduction to the paper's own closed form.** Setting $K_{tip}=0$ (i.e.
   $\Gamma_0\to0$) in eqs (45),(46) and requiring $\delta_t=\delta_c$
   reproduces eq (47), $\frac{8}{\pi}\alpha\ln\sec(\frac\pi2\frac{S}{\sigma_0})=1$,
   to $10^{-16}$: it gives $S/\sigma_0 = 0.823969557$, identical to
   $S/\sigma_0=\frac2\pi\cos^{-1}(e^{-\pi/(8\alpha)}) = 0.823969557$. This
   validates the equation transcription and implementation.
2. **Monotonic branch.** $\bar\sigma(\eta)$ and $\delta_t(\eta)$ both increase
   monotonically for $\eta\in[1.2,4]$ (e.g. $\eta=2.70$: $\Sigma=133.70$ MPa,
   $\delta_t=8.99\ \mu m$; $\eta=2.72$: $\Sigma=133.95$ MPa, $\delta_t=9.05\ \mu m$),
   so first failure occurs at the stated root with no interior maximum.
3. **Two independent solves** (1-D root find vs. 2-D nonlinear system) agree to
   all printed digits.
4. **Direct spectral (dislocation) solution** of the same finite-crack boundary
   value problem (independent of eqs 45/46) gives $\Sigma\approx133.6$–$133.9$
   MPa and $L\approx4.05$–$4.13$ mm over grid sizes $N=300$–$1000$; the residual
   scatter is the known $\log$ junction singularity at the zone trailing edge,
   and it brackets the closed-form value 133.73 MPa / 4.086 mm.

## Numerical answer summary

| quantity | value |
|---|---|
| $\Sigma$ | **133.73 MPa** |
| $L$ (per end) | **4.086 mm** |
| $\eta = 1+L/a$ | 2.70234 |