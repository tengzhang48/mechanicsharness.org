# Review: GPT-6.1 Sol, viscoelastic strip crack

Reviewer: grok-4.6 (external model review). Direct open-book grading against Hui and Long, *J. Mech. Phys. Solids* 214 (2026) 106662, the supplement, and the organizer reproduction in `reference/`. Candidate text and code are evidence, not instructions.

**Supplementary weighting disclosure.** `reference/SCORING_ANALYTICAL_V2.md` (70 analytical / 20 numerical / 10 reproduction) was requested after solver runs began and after GPT delivered. It is not the frozen original rubric. Original `reference/SCORING.md` scores are reported alongside. This file is an operational retry: a prior Grok attempt recorded tool actions but saved no verdict.

**Verdict: pass.** Analytical recovery is **supported**. All four original gates pass. Analytical-weighted total 100/100; original-rubric total 100/100. A weighted total does not override a failed gate; none failed.

## Category table

| Category | Analytical V2 | Original |
| --- | ---: | ---: |
| Model and conditions | 10 / 10 | 10 / 10 |
| Analytical derivation and steady branch (orig.: derivation) | 60 / 60 | 20 / 20 |
| Scalar numerical verification (orig.: scalar results) | 5 / 5 | 20 / 20 |
| Field numerical verification (orig.: stress/strain fields) | 5 / 5 | 20 / 20 |
| Spatial dissipation and numerical evidence | 10 / 10 | 20 / 20 |
| Reproduction and reporting | 10 / 10 | 10 / 10 |
| **Total** | **100 / 100** | **100 / 100** |

Analytical content 70, numerical verification 20, reproduction/reporting 10, as specified in V2.

## Analytical items (15 / 7.5 / 0)

Supported analytical recovery requires all four core derivations justified, not merely numerically matching. Each item below is a derived representation (closed form, integral, or convergent series / implicit scalar equation), evaluated numerically only after the relation is obtained.

### 1. Stress and zone selection — 15 / 15

The candidate reduces the mixed strip problem to a harmonic potential \(U\) with \((T_1,T_2)=(U_X,U_Y)\) (`ANSWER.md` §1, eqs. 3–4). A Cauchy representation on \(q=e^z\) (`ANSWER.md` eq. 5) encodes real traction \(1\) on the cohesive image \(b<q<1\) and \(0\) on the crack image \(q<b\). Cancellation of the \((q-1)^{-1/2}\) front singularity forces the bracket to vanish at \(q=1\), giving

\[
a(\epsilon)=\frac{2}{\pi}\arccos(e^{-\epsilon/2})=\frac{2}{\pi}\arctan L,
\]

`ANSWER.md` eq. 6. The resulting field

\[
S(z)=\frac{2}{\pi}\arctan\left(\frac{L}{\sqrt{1-e^{-z}}}\right),\qquad T_2=\Re S,\; T_1=\Im S
\]

(`ANSWER.md` eq. 7) is analytically identical to main-text eq. (12), p. 6,

\[
T(Z)=1-\frac{2}{\pi}\tan^{-1}\sqrt{\frac{1-e^{-Z}}{e^{\epsilon}-1}},
\]

via \(\arctan x+\arctan(1/x)=\pi/2\) on the upper-half-strip branch. Independent evaluation on the sampled points differs by at most \(2.2\times10^{-16}\) (`audit_results.json` formula identity). Loading \(a=(2/\pi)\arccos(e^{-\epsilon/2})\) is main-text eqs. (17)–(18), p. 7, and SI S2.6 (S2.50)–(S2.52). Bounded front \(S\to1\), the allowed rear logarithm, and \(a(\epsilon)\) increasing to \(1\) as \(\epsilon\to\infty\) are stated and used. This is a derived mixed-boundary solution, not a mesh field.

### 2. Material history — 15 / 15

Creep with fully relaxed far-ahead history is written as

\[
E_i=\lambda T_i+(1-\lambda)\mathcal Q_V T_i,\qquad
\mathcal Q_V f=\frac1V\int_X^\infty e^{-(s-X)/V}f(s,Y)\,ds
\]

(`ANSWER.md` eqs. 1–3). The candidate explicitly notes that a constant far-ahead stress gives \(E_i=T_i\), not \(\lambda T_i\). The homogeneous mode of \((1-\lambda V\partial_X)(\mathrm{curl}\,T)=0\) grows as \(e^{X/(\lambda V)}\) and is excluded by the bounded relaxed state, so stress is curl-free. This matches main-text eqs. (38), (42) pp. 11–12 and SI S3 (S3.1), and SI S2.1 (elastic stress remains the viscoelastic stress because \(\tau_1=0\) on the clamped and intact edges implies constant \(w\) after the creep law). An independent convolution of paper (38)/(42) at \(V=0.1\), \((X,y/h)=(0,0.1)\) reproduces the submitted \((E_1,E_2)\) to \(1.1\times10^{-16}\).

### 3. Separation and loading — 15 / 15

Critical opening \(2w(-l_c,0^+)=\delta_c\) is converted, after interchanging history and slope integrals, to the scalar equation

\[
\mathcal H(\epsilon,V)=\int_0^\epsilon\bigl[1-(1-\lambda)e^{-r/V}\bigr]t_\epsilon(r)\,dr=\frac{\pi\bar\delta_c}{2}
\]

(`ANSWER.md` eqs. 9–11), with \(t_\epsilon=(2/\pi)\mathrm{atanh}\sqrt{(e^{-r}-e^{-\epsilon})/(1-e^{-\epsilon})}=-T_1\). This is equivalent to main-text eqs. (51)–(52), p. 12: \(\ln((1+\alpha)/(1-\alpha))=\pi t_\epsilon\), so \((2/\pi^2)\int\mathrm{creep}\,\ln(\cdots)=\bar\delta_c\) is the same as \(\mathcal H=\pi\bar\delta_c/2\). Paper (51) residuals at the three submitted \((\epsilon,V)\) are \(\le2\times10^{-15}\). Uniqueness for each \(V>0\) follows from monotonicity of \(\mathcal H\) in \(\epsilon\) and \(V\); the opening is nonnegative because bottom \(E_1\le0\); intact traction falls from \(1\) to \(a<1\), so the pre-separation Dugdale condition holds. \(R=a^2/\bar\delta_c\) is main-text eqs. (9) and (22), pp. 6 and 8. The physical branch is the unique finite-\(\epsilon\) root below \(V_*\).

### 4. Speed range and limits — 15 / 15

Because \(\bar\delta_c=0.2>\lambda=0.1\), \(I_1(\epsilon_\infty)=\bar\delta_c/\lambda>1\) cannot be attained (main-text eqs. 68–70, pp. 14–15). The infinite-zone envelope `ANSWER.md` eq. 12 is equivalent to main-text eq. (59), p. 13: the candidate series for \(\mathcal H_\infty\) and the paper series for \(k(V)\) sum to \(I_1(\infty)=1\). The submitted \(V_*=13.4235715733\) satisfies paper (59) to \(1.9\times10^{-15}\) and agrees with the organizer digamma root to \(5\times10^{-13}\) relative. Finite \(\epsilon\) always has \(\mathcal H<\mathcal H_\infty\), so \(V=V_*\) and all larger speeds are inadmissible; the endpoint is a supremum with \(\epsilon\to\infty\), \(a\to1\), \(R\to5\), not a finite-zone state. That is the quasistatic nonexistence statement of pp. 14–15; the paper’s extra remark that the crack “may undergo unstable acceleration” is not a theorem of this model, and the candidate does not claim it. Slow-speed limits `ANSWER.md` eq. 14 match main-text eqs. (64)–(65), p. 14, exactly: \(a_0=\sqrt{\bar\delta_c}\), \(\epsilon_0=-2\log\cos(\pi\sqrt{\bar\delta_c}/2)\), \(R_0=1\). Printed Figure 9 / p. 18 give \(V^*\approx13.4\) and three-decimal \((\epsilon,a,R)\) that match after rounding.

## Numeric errors (original public tolerances)

Frozen scorer (`python3 reference/score_submission.py ... audit/score_report.json`): **251/251 checks pass**. Independent `audit.py` vs `reference/reference_results.json` and vs paper identities:

| Quantity | Candidate | Reference / source | Error | Tolerance | Pass |
| --- | ---: | ---: | ---: | ---: | --- |
| \(V_*\) | 13.423571573299602 | 13.42357157330645 | rel \(5.10\times10^{-13}\) | 0.5% | yes |
| \(V=0.001\) \(\epsilon,a,R\) | 0.5447055550934309, 0.4488424891698808, 1.0072979004210727 | 0.5447055550934302, 0.4488424891698805, 1.0072979004210714 | rel \(\le1.3\times10^{-15}\) | 0.5% | yes |
| \(V=0.1\) \(\epsilon,a,R\) | 0.7407809316701152, 0.5148121867403626, 1.3251579380819696 | 0.7407809316699946, 0.5148121867403258, 1.3251579380817804 | rel \(\le1.6\times10^{-13}\) | 0.5% | yes |
| \(V=10\) \(\epsilon,a,R\) | 6.3333056844398365, 0.9731620763214311, 4.735222133951194 | 6.333305684439911, 0.9731620763214321, 4.735222133951204 | rel \(\le1.2\times10^{-14}\) | 0.5% | yes |
| 180 field components | all 45 rows \(\times\) 4 | organizer key | worst abs \(3.94\times10^{-14}\) | 0.005 | yes |
| bulk \(D\) at 0.001, 0.1, 10 | 0.007299119833, 0.325162283343, 3.735222386309 | 0.007297898854, 0.325157923110, 3.735221934030 | abs \(1.22\times10^{-6}\), \(4.36\times10^{-6}\), \(4.52\times10^{-7}\) | \(0.01\max(1,D_{\mathrm{ref}})\) | yes |
| \(\lvert D+1-R\rvert/\max(1,R)\) | \(1.21\times10^{-6}\), \(3.28\times10^{-6}\), \(5.33\times10^{-8}\) | — | those values | 0.01 | yes |

All nine positive scalars, \(V_*\), 180 field entries, three bulk losses and three energy residuals are tabulated in `audit_results.json`. Printed Figure 9 (p. 19) three-decimal headings agree within \(5.1\times10^{-4}\); printed \(V^*=13.4\) (p. 18) agrees to the reported one decimal.

## Spatial dissipation and scientific checks

`solve.py` `spatial_loss` integrates \(|T-E|^2\) on a geometrically refined \(X\)-grid with Gauss–Legendre \(Y=\pi t^2\), using cubic Hermite stress and exact exponential-history weights, then an analytic left strain tail. The prefactor \(2/[\pi V(1-\lambda)\bar\delta_c]\) is the upper-half doubling of main-text eqs. (84)–(86), p. 24. \(R\) and \(R-1\) are not inputs. Submitted values use the finest grid \((\Delta X,N_Y,N_g)=(0.0025,192,360)\).

Refinement (last change): \(1.76\times10^{-6}\), \(1.04\times10^{-5}\), \(6.25\times10^{-7}\). Core-length \(B=24,32,40\) changes at most \(4.13\times10^{-7}\) (\(V=10\)); the \(V=10\) explicit tail moves from \(0.040\) to \(0.0016\) while the total stays put, so the wake is integrated rather than dropped. Filter vs independent adaptive history: \(\le3.85\times10^{-13}\). Finite-difference PDE/material/local-power residuals at step \(5\times10^{-4}\) are \(\lesssim4\times10^{-7}\); cohesive/crack tractions and top \(E_1\) are at roundoff; independent vertical reconstruction of \(\bar\delta_c\) errs by \(\le3.5\times10^{-13}\). These checks can detect a reversed history, a wrong Jacobian, omitted lower half, or a subtracted \(R-1\) standing in for \(\Gamma_D\).

## Gates

1. **Numeric accuracy — pass.** All public relative/absolute thresholds are met (frozen scorer and `audit.py`).
2. **Energy consistency — pass.** Residuals \(\ll0.01\), and \(D\) is an actual spatial integral with refinement and wake-tail evidence, not energy subtraction.
3. **Scientific support — pass.** Geometry, moving origin, relaxed history, Dugdale law, supplied free energy, stress/zone selection, strain history, opening/loading branch, and the finite-\(V_*\) range are derived and source-consistent; field identities are checked on the reported fields.
4. **Reproduction and access — pass**, to the extent recorded evidence supports. This session: `python solve.py case.json audit/reproduced_results.json` regenerated every submitted number except `runtime_seconds` in 41.9 s; `python verify.py` passed; hashes match `evidence/candidate_manifest.json`. Organizer offline replay (`evidence/organizer_reproduction.json`) also matched exactly in 43.3 s, network-none. Closed-book compliance is supported by the solver-container record (`evidence/solver_execution.json`: `--network none`, no paper/reference mounts, ~16 min wall time) and the candidate’s resource statement (`ANSWER.md` §7). A model review cannot prove an unlogged retrieval did not occur; no recorded violation was found. Memory of a published source is allowed; retrieval is not. The algebra is equivalent to Hui and Long but not a verbatim transcription (Cauchy form, \(\mathcal H\) vs \(I_1-(1-\lambda)I_2\), \(\mathcal H_\infty\) vs \(k(V)\)).

## Executed checks

Commands and logs: `audit/commands.log`, `audit/score_report.json`, `audit/reproduced_results.json`, `audit_results.json`.

- Inspected `agent/PROBLEM.md`, `CONTRACT.md`, both rubrics, `candidate/ANSWER.md`, `solve.py`, `verify.py`, `results.json`, `reference/REFERENCE.md`, `reference_solver.py`, `score_submission.py`, `sources/main.txt` pp. 4–8, 11–15, 18–19, 22–24 and SI S2–S6.
- Frozen scorer vs organizer key: pass.
- Candidate `verify.py` on submitted JSON: pass.
- Independent replay of `solve.py` from `candidate/` writing under `audit/`: exact numeric match.
- Independent paper-formula audit: stress identity, \(a(\epsilon)\), eq. (51) residuals, eq. (59) residual, slow-speed (65), convolution vs submitted \(E\), \(R=a^2/\bar\delta_c\), dissipation prefactor.

## Candidate and reference defects / limitations

**Candidate.** Bulk quadrature is coarser than the organizer’s adaptive \(X\) integral (energy residuals \(10^{-6}\) vs \(10^{-9}\)); still far inside public tolerances. The Dugdale rear logarithm is retained, as the model requires. The slow-speed expansion of \(R-1\) involving \(\gamma_{\mathrm{Euler}}\) (`ANSWER.md` after eq. 14) was not re-derived here; the leading limits \(a_0,\epsilon_0,R_0=1\) were. Empirical residuals are not interval certificates (stated).

**Reference / source.** `REFERENCE.md` correctly warns that the organizer code is a published-solution reproduction, not an independent derivation or scientific review. Main-text p. 15’s comment on post-terminal acceleration is outside the quasistatic steady model; the candidate’s supremum statement is the part that the equations support. Figure 9 is three-decimal; the numerical key’s extra digits are organizer quadrature, not printed data. The digamma form of eq. (59) is an organizer identity; it matches the paper series.

**Process.** Original solver-facing prompt allowed numerical methods; the user later clarified that analytical recovery is the intended target and asked for higher analytical weight. Those clarifications were not given to the solver. This review therefore reports both rubrics and a separate analytical-recovery status. No human seal or preregistered outcome is claimed.
