# Independent check of the poroelastic small-time tip volume

**Inference / recommendation:** retain the tip term **1.086165**. A direct finite-crack solve, using no factorization
values, supports **1.09 ± 0.01**. A newly implemented numerical Cauchy solve gives **1.0861646**, with an estimated
numerical uncertainty of **±0.00002** within the shared tip model. Thus the independent physical calculation confirms the
effect and its size, not six decimal places. These are sensitivity estimates, not certified bounds. No key was edited.

## Evidence and source equations

**Read:** start at `../FIXES.md:9–43`; the frozen task is `../task/PROBLEM.md:3–43`; the current key states the contested
decomposition at `../key/KEY.md:17–20,75–84`. The external audit is
`audit_2026-10-25/AUDIT.md:38–53` (the table and its qualifications).
The external `FE_CHECK.md`, at
`key/FE_CHECK.md`,
lines 70–131, especially “FE with the face layer resolved”, reports A ≈ 2.2–2.4; it does not resolve the digits.

**Read / visually checked:** source file
`<source paper PDF, not included>`,
Detournay and Cheng (1991), PDF pages 5–7, printed pp. 1649–1651, Eqs. (15)–(31). Equations below are transcribed from
those rendered pages, with subscripts/typesetting regularized; they are not quotations from the key:

\[
 V_c^{(1)}(0^+)=\frac{\pi L^2(1-\nu_u)}G,\qquad
 \lim_{t\to0^+}p=-\frac{\nu_u-\nu}{\alpha(1-2\nu)}
 \lim_{t\to0^+}(\sigma_{xx}+\sigma_{yy}) .\tag{15,16}
\]
\[
 \lim_{t\to0^+}p=\frac{2(\nu_u-\nu)}{\alpha(1-2\nu)}
 \quad(y=0,\ -L\leq x\leq L),\tag{17}
\]
\[
 V_c^{(1)}(t)=\frac{\pi L^2(1-\nu_u)}G+
 \frac{16L^2(1-\nu_u)(\nu_u-\nu)}{G\sqrt\pi(1-\nu)}
 \sqrt{ct/L^2},\qquad \sqrt{ct/L^2}\ll1 .\tag{18}
\]
**Read:** the text immediately above (18) explicitly approximates the crack-surface diffusion as one-dimensional.
PDF p. 6 / printed p. 1650 calls its SIF reasoning a “heuristic argument”; PDF p. 7 / printed p. 1651 displays
\(K_I^{(1)}(0^+)=[(1-\nu_u)/(1-\nu)]\sqrt{\pi L}\) in (29). The same page gives
\(V_c^{(2)}(t)=V_l^{(1)}(t)\) in (31), following the reciprocity integral (30).

**Inference:** (18) supplies the trusted flat-face term. It is a simplified solution, not a paper error. The original
benchmark key's treatment of that term as the complete coefficient was our error. No claim about the authors' intent
beyond their printed approximation and heuristic is needed.

## 1. Re-derive and test the boundary operator

**Computed / derivation checked by `airy_mode_check.py`:** set G=L=c=1, d=νu−ν, a=α(1−2ν), b=a/(1−ν),
η=d/(1−νu)=1/3, β=η/b and C=G/[2(1−ν)]. Compatibility and fluid mass balance give, with T=σxx+σyy,
\(h=T+bp\),
\[
 \nabla^2h=0,\qquad \nabla^2p=q^2(p+\beta h),\qquad q=\sqrt{s/c}.
\]
For Fourier wavenumber k=|ξ| and m=√(k²+q²), use decaying Airy modes
\(\Phi=(A+By)e^{-ky}+De^{-my}\),
\(p=Pe^{-my}-\beta He^{-ky}\). Enforcing shear-free symmetry and integrating the supplied yy strain gives
\[
 S=-Ck\left[1+\eta\frac{m-k}{m+k}\right]w-\frac{bk}{m+k}p_0,
 \qquad Q=-mp_0+2\beta Ck(m-k)w .
\]
Here w is full opening, S=σyy, Q=∂yp and p0=p(x,0). The script checks compatibility, zero shear, h,
both mass-balance modes, opening and both surface relations symbolically; **all eight residuals are exactly zero**.

**Computed / implementation:** eliminate p0 and use F=−Q/√λ, λ=2βCq²/b. The resulting symmetric operator is
\[
 \binom{-S}{p_0/\sqrt\lambda}=
 \begin{pmatrix}
 Ck[1+\eta-2\eta k^2/(m(m+k))]&\sqrt{2\eta C}\,qk/[m(m+k)]\\
 \sqrt{2\eta C}\,qk/[m(m+k)]&1/m
 \end{pmatrix}\binom wF .
\]
On |x|<1 impose −S=1, p0=0; outside impose w=Q=0. In particular, pressure is **not** set to zero ahead of the tip.
This solves the full finite slit, including its face diffusion and both tips.

## 2. Route A: direct finite-crack solve, without factorization

**Computed by `finite_crack_check.py`:** apply this operator by FFT on a periodic interval of length W and solve its
restriction to the slit with preconditioned conjugate gradients. A same-grid undrained elastic solution is subtracted
before extracting the resolvent volume slope. For U(s)=s V̂(s),
\[
 A(q)=\frac{2q}{\sqrt\pi}[U(s)-U(\infty)] .
\]
The factor 2/√π follows from the transform of √t. Tip intensity is obtained from the opening ratio to that same-grid
elastic solution, multiplied by (1−νu)/(1−ν)=0.75; the mesh is refined before the large-q extrapolation.

**Computed / calibrations first:** `face_calibration()` numerically integrates erfc(y/2) on the half-line,
then applies the pressure/reciprocity multiplier. It prints **1.3540550005**, matching 1.354055.
The elastic volume tends to 0.6π as W increases (the periodic image bias remains at fixed W).
The intensity extrapolations below independently recover **0.80743–0.80754**, within **0.0002** of trusted 0.807593.
No factor entries, energy identity or target A enter the direct solver or either extrapolation.

**Computed / refinement:** W=16 uses N=262144,524288,1048576; W=32 and W=64 each use
N=524288,1048576,2097152. q=80,160,320,640,1280. For W=32 the finest q dx is at most 0.01953125;
for W=64 it is at most 0.0390625. All reported direct residuals are below 2×10⁻¹¹.
First-order extrapolation from the last two meshes is followed by W⁻² extrapolation from W=32,64.
Including W⁻⁴ using all three boxes changes A(q) by at most 3.3×10⁻⁶.

**Computed by `summarize_direct.py`, infinite-box, zero-mesh estimates:**

| q | A(q) | calibrated K(q) |
|---:|---:|---:|
| 80 | 2.1153300 | 0.8276050 |
| 160 | 2.2054834 | 0.8216896 |
| 320 | 2.2715176 | 0.8175307 |
| 640 | 2.3194641 | 0.8146001 |
| 1280 | 2.3540424 | 0.8125208 |

**Computed / extrapolation sensitivity:** polynomials of degrees 1,2,3 in q⁻¹ᐟ², using all five or the last four
points, give A∞=2.431965–2.440006. Quadratic/cubic fits give 2.439375–2.440006 and
tip=1.085320–1.085951. Even the deliberately lower-order linear fits are retained in the uncertainty estimate.

**Inference:** report this route conservatively as A=2.44±0.01, tip=1.09±0.01. Finite-q extrapolation dominates;
the small algebraic residual is not a continuum error bound. The pilot range q=20–160 did not reach the asymptote:
an initial q⁻¹/log(q)/q fit was rejected after observing q⁻¹ᐟ² approach and testing the trusted K. It is not evidence
against the key. This is a Laplace-domain boundary calculation, not a new time-stepping FE run.

## 3. Route B: numerical Cauchy factor and its limits of independence

**Computed by `cauchy_check.py`:** remove the known scalar √(±iξ) and √(1±iξ) factors from the surface symbol.
For the remaining Hermitian matrix N(ξ), solve the linear Riemann–Hilbert equation
\(X+P_-[(N-I)X]=I\), with X lower analytic and X(∞)=I. P− is a numerical Cauchy/Riesz projection on
ξ=scale tan(θ/2); GMRES solves for X−I. The lower factor is X⁻¹. This code imports no archived solver and uses no
saved factor values. It is distinct from the archived nonlinear Wilson calculation.

**Read / qualification:** the archived `../runs/MP-37/workspace/factor_check.py:24–44` also solves a projected inverse-factor
equation. Thus this is a new numerical implementation in the same factorization family, **not independent physical
validation**. Its symbol and the matching/volume identities are shared with the archived analysis in
`../grading/review_B_K10-codex/output/VERDICT.md:42–332`, sections 1–3. Route A supplies the independent check requested here.

**Computed / matching:** with n11,n21 the first column of the normalized lower factor at zero,
\[
 K_0=\frac{n_{11}}{1+\eta},\quad
 A_{face}=\frac{16d}{(1+\eta)\sqrt\pi},\quad
 A_{tip}=\frac{4d}{(1+\eta)\sqrt\pi}\pi\sqrt{2/\eta}\,n_{21}.
\]
These identities are tested against the directly computed volume above; that comparison does not assume their answer.

**Computed**, N=4096,16384,65536, extrapolating the zero-frequency square-root cusp with N⁻¹ᐟ² and N⁻¹:

| map scale | fit | K0 (calibration) | tip | faces + tips |
|---:|---|---:|---:|---:|
| 0.01 | linear | 0.807592984 | 1.086164543 | 2.440219544 |
| 0.01 | quadratic | 0.807592993 | 1.086164578 | 2.440219579 |
| 0.1 | linear | 0.807592759 | 1.086164237 | 2.440219238 |
| 0.1 | quadratic | 0.807592967 | 1.086164553 | 2.440219553 |

**Computed / additional diagnostic:** the configurational identity n11²+n21²=1+η, with **trusted K0 supplied as an
input**, yields tip=1.086164605. This is conditional consistency, not a third independent measurement of the tip term.

**Inference:** the Cauchy numerical uncertainty ±0.00002 comfortably exceeds the observed extrapolation spread, but
does not bound model or matching errors. Route A independently limits those errors at roughly the percent level.
The six-digit rounding is consistent; a six-digit physically independent certification remains unavailable.

## Reproduction and disposition

**Computed / commands**, from this folder, using system Python with numpy, scipy, sympy, mpmath:

```sh
python3 airy_mode_check.py
python3 cauchy_check.py
python3 finite_crack_check.py
python3 finite_crack_check.py --sizes 262144 524288 1048576 --box 16 --qs 80 160 320 640 1280
python3 finite_crack_check.py --sizes 524288 1048576 2097152 --box 64
python3 summarize_direct.py
```

**Read / reproducibility:** the three finite-crack commands produced `finite_crack_results.jsonl`,
`finite_crack_box16.jsonl`, `finite_crack_box64.jsonl` respectively, by stdout redirection.
`cauchy_results.jsonl` and `direct_summary.jsonl` preserve the other numeric tables. Scripts print their results and
do not modify the task, key, checker or saved runs; the summary script reads these small checked-in diagnostics.

**Recommendation for Teng:** close the outstanding *independent size check* at the stated uncertainty; keep 1.086165
as the numerical factorization reference. Keep the qualification that its last digits rely on the common analytical
matching model. Any change to key wording, grading tolerances or the public case page is a separate approval decision.
