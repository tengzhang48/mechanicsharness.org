# Grading: submission C

Fresh review of anonymous submission `C/` against `PROBLEM.md` and the corrected key in `REFERENCE.md`. No earlier model verdicts were used. Supplied files under `C/` were not modified. Independent checks used copies under `scratch/`.

**Overall: correct.**

The nine reported numbers all lie inside the stated tolerances (1% relative where nonzero; absolute 0.01 for `K2_0`). The physical boundary-value problem is the one asked for: total (not effective) face traction, prescribed face pore pressure, undrained outer field with a drained inner tip, and a genuinely coupled loading-2 problem whose long-time and far-field limits are nonuniform in two dimensions. Closed-form items are derived, not merely fitted. The two quantities that have no elementary closed form (`K1_0`, the tip part of `V1_short_coeff`) are obtained from a documented Laplace-domain mixed BVP, with saved residuals, convergence tests, and an extrapolation whose structure matches the inner-layer scaling. Numerical agreement with the key is not treated as a proof; the formulation, identities, and representative reruns are.

The original wording of “just after loading (\(t=0^+\))” leaves the order of the tip and time limits implicit. That ambiguity is preserved below: an undrained-outer reading of 1 is a legitimate reading of the sentence, and is not marked wrong. The submission reports the drained-tip limit as `K1_0` and records the other reading explicitly, which is the handling the corrected key asks for.

---

## 1. What was asked

Plane-strain Biot crack, faces \(|x|\le L\), \(y=0\).

- Loading 1: \(\sigma_n=-H(t)\), \(p=0\) (total normal traction and drained faces).
- Loading 2: \(\sigma_n=0\), \(p=H(t)\).

Report \(K\) in units of \(\sqrt{\pi L}\) and \(V\) in units of \(L^2/G\).

Corrected key (`REFERENCE.md`, log P374), \(\nu=0.2\), \(\nu_u=0.4\), \(\alpha=0.89\):

| Quantity | Corrected key | Paper (on record, not the key) |
|---|---:|---:|
| `K1_0` | **0.807593** | 0.75 (heuristic, Eq. 29) |
| `K1_inf` | 1 | 1 |
| `V1_0` | \(\pi(1-\nu_u)=1.884956\) | same |
| `V1_inf` | \(\pi(1-\nu)=2.513274\) | same |
| `V1_short_coeff` | **2.440220** = 1.354055 (faces) + 1.086165 (tips) | 1.354 (faces only, Eq. 18) |
| `p_faces_0` | \(2(\nu_u-\nu)/[\alpha(1-2\nu)]=0.749064\) | same |
| `K2_0` | 0 | 0 |
| `K2_inf` | \(-\eta=-0.33375\) | same |
| `V2_inf` | \(-\pi\alpha(1-2\nu)/2=-0.838805\) | same |

\(\eta=\alpha(1-2\nu)/[2(1-\nu)]\). Tolerance: 1% relative if nonzero; `K2_0` within 0.01.

The paper’s `K1_0=0.75` is neither limit of the inner problem. The paper’s short-time volume omits the tip regions. Those two corrections are the decisive points of the problem.

---

## 2. Physical conditions (the critical ones)

### 2.1 Total traction and pore-pressure face conditions — correct

The constitutive law is given with tension-positive total stress \(\sigma_{ij}\) and pore pressure \(p\). The face data are \(\sigma_n=\sigma_{yy}\) and \(p\), not an effective-stress combination \(\sigma_{yy}-\alpha p\).

The solver imposes exactly those mixed data in the Laplace domain (`C/code/crack_solver2.py` 104–113):

```python
for name, (sc, pc) in {'L1': (-1.0/s, 0.0), 'L2': (0.0, 1.0/s)}.items():
    rhs[0] = sc*np.pi/2; rhs[N] = pc*np.pi
```

with \(|x|<L\): \(\sigma_{yy}\) and \(p\) prescribed, \(|x|>L\): opening \(U=0\) and flux \(Q=0\) (`ANSWER.md` 250–252). The half-plane kernel is built from \(\sigma_{xy}=0\), \(u_y=U\), \(q_y=Q\) on \(y=0\), returning \(\sigma_{yy}\) and \(p\) (`symbols.py` 58–65; remainders in `crack_solver2.py` 39–44).

Loading 2 is therefore not replaced by a naive Biot effective tension \(-\alpha\) on the faces. The long-time coefficient is \(\eta=0.33375\), not \(\alpha=0.89\). That is the coupled back-stress, not a constitutive rewrite of the BC.

### 2.2 Undrained outer vs drained inner tip; order of limits — correct, ambiguity preserved

The problem statement says only “just after loading (\(t=0^+\))”. Two limits exist (`REFERENCE.md` 24–31):

1. **Outer / simultaneous undrained snapshot** (\(t=0\) with \(\zeta=0\) everywhere, then the elastic \(K\)-field): \(K=\sigma\sqrt{\pi L}=1\) in the requested units. The opening is the undrained Sneddon ellipse, which is also \(V(0^+)\).
2. **Inner tip at any \(t>0\)**, then \(t\to 0^+\) (\(r\to 0\) at fixed \(t>0\)): a self-similar drained zone of size \(\sqrt{ct}\) sits at each tip, and \(K_{\mathrm{tip}}=a(\gamma)\,K_\infty\) with \(\gamma=(\nu_u-\nu)/(1-\nu)\). For the given data, \(a(0.25)=0.80759\). This is the corrected key.

Neither limit is 0.75. The paper’s heuristic (undrained opening inserted into the drained \(K\)–opening relation, Eq. 29) gives \(\sqrt{(1-\nu_u)/(1-\nu)}=0.75\), which is a lower bound from a Laplace-domain \(J\)-integral only after the fluid tip term is dropped (`ANSWER.md` 196–204). The submission treats 0.866 as an *upper* bound on \(a\) and 0.75 as neither bound nor limit.

The submission states both readings (`ANSWER.md` 35–39, 337–340; `results.json` `K1_0` and `supplementary.K1_0_instantaneous_undrained_t_exactly_0`). It reports **0.80759** as `K1_0` because the tip is drained for every \(t>0\), and records 1 as the instantaneous undrained outer value. That is the right way to preserve the implicit order of limits: 1 is not marked as an error here, and 0.80759 is the value the corrected key uses.

The inner/outer split is implemented, not just described:

- far field: undrained \(K\)-field with \(K_\infty=\sigma\sqrt{\pi L}\);
- faces: traction-free and \(p=0\);
- \(K\) extracted with the **drained** conversion \(K/\sqrt{\pi L}=G\sum(n+1)a_n/(1-\nu)\) (`ANSWER.md` 265; `crack_solver2.py` 113);
- the combination \(W=s(\tilde K_1+2\beta\tilde K_2)\) cancels the uniform face-pressure mismatch \(2\beta\), leaving a remainder that approaches \(a(\gamma)\) as \(s^{-1/2}\) (`final_extrapolation.py` 23–32; `ANSWER.md` 211–213).

Four parameter sets with the same \(\gamma=0.25\) are claimed to give identical \(W\) (`ANSWER.md` 311–312; `check_gamma_scaling.py`). That is the right invariance.

### 2.3 Loading 2: coupled problem, nonuniform long-time and far-field limits in 2-D — correct

A uniform pore pressure in an infinite 2-D body, stress-free at infinity, is *not* a uniform infinite-plane state with \(p\equiv 1\). Stresses and \(p\) vanish at infinity for all \(t\) (problem statement). In two dimensions the pressurized region grows as \(\sqrt{ct}\) and the deficit near the crack decays only as \(1/\ln t\).

The submission’s argument (`ANSWER.md` 97–121):

- a large disc of uniform \(p\) in an infinite plane has uniform interior stress \(\sigma_{xx}=\sigma_{yy}=-\eta p\), independent of radius (thermoelastic analogy);
- uniform \(p\) alone is free swelling: no stress, no opening;
- the crack faces must remain free of *total* normal traction, so the correction is a Griffith crack with face tension \(+\eta\);
- therefore \(K_2(\infty)=-\eta\sqrt{\pi L}\) and \(V_2(\infty)=-\pi(1-\nu)\eta L^2/G=-\pi\alpha(1-2\nu)L^2/(2G)\).

That is the coupled long-time problem, not “put \(p=1\) everywhere and use undrained or drained moduli.” The far-field remains unloaded; the \(-\eta\) plateau is a slowly expanding near-crack field. The approach is logarithmic, with capacity constant of a segment of length \(2L\):

\[
s\tilde K_2=-\eta\bigl[1-1/(\ln(1/s)+4\ln 2-2\gamma_E)\bigr],
\]

and in time \(K_2\approx-\eta[1-1/(\ln(16\tau)-\gamma_E)]\). The small-\(s\) fit (`small_s.log`) gives \(X_\infty=-0.333750\), \(D=1.6181\), and \(4\ln 2-2\gamma_E=1.618157\) independently. At \(\tau=10^4\), saved Talbot data still have \(K_2/K_2(\infty)\approx 0.914\). Treating a moderate-time numerical value as \(K_2(\infty)\) would have been wrong; they did not.

\(K_2(0^+)=0\) follows from the inner problem being driven only by a layer of thickness \(\sqrt{ct}\), so \(K_2=O(t^{1/4})\to 0\). Sign: swelling walls close the crack, \(K_2<0\) for \(t>0\).

---

## 3. Part (a) — correct

### \(V_1(0^+)\) and \(V_1(\infty)\) — correct (closed form)

Undrained then drained Griffith volumes:

\[
V_1(0^+)=\pi(1-\nu_u)\sigma L^2/G,\qquad V_1(\infty)=\pi(1-\nu)\sigma L^2/G.
\]

Reported 1.884956 and 2.513274. Ratio on the agent’s own numbers: \(V_1(\infty)/V_1(0^+)=1.333333\), equal to \((1-\nu)/(1-\nu_u)\) (principle check in `REFERENCE.md` 105–109).

### \(K_1(\infty)\) — correct (closed form)

Drained Griffith value \(\sigma\sqrt{\pi L}\), independent of material constants. Reported 1.

### \(K_1(0^+)\) — correct as the drained-tip limit; the other reading is recorded

Reported 0.80759 against key 0.807593 (relative error \(4\times 10^{-6}\)).

Reasoning is the inner self-similar problem, not the paper heuristic. Weak-coupling expansion \(a=1-3\gamma/4+O(\gamma^2)\) is derived from the erf regularisation of the tip suction and a Laplace-domain \(J\)-integral that keeps the fluid term (`ANSWER.md` 196–209). At \(\gamma=0.01\), saved data give \((1-a)/\gamma=0.75069\to 3/4\). No closed factorization of the \(2\times 2\) Wiener–Hopf kernel is claimed.

**This number is a numerical BVP result, not an analytic proof.** Agreement with the key, with the Yu/Atkinson–Craster numerical value quoted in the key, and with the submission’s own \(J\)-bound, does not replace a factorization. What *is* shown is that the quantity being extrapolated is the drained-tip factor of the mixed problem, that the face layer has been subtracted, and that the large-\(s\) remainder behaves as \(s^{-1/2}\).

### Short-time volume — correct (faces exact + tips numerical)

The paper’s Eq. 18 is the face-layer piece only. The submission splits \(A=A_{\mathrm{face}}+A_{\mathrm{tip}}\) at the same order \(O(\sqrt{t})\).

**Faces (exact).** Induced undrained pressure next to the faces \(p=2\beta\sigma\), incompatible with \(p=0\) on the crack, dissipated by 1-D diffusion \(\Delta p=-2\beta\sigma\,\mathrm{erfc}(y/2\sqrt{ct})\). Reciprocal identity (R) (`ANSWER.md` 71–80) converts \(\zeta\) in the layers into volume. Result:

\[
A_{\mathrm{face}}=\frac{16}{\sqrt{\pi}}\gamma(1-\nu_u)=1.354055,
\]

which is the paper’s coefficient. A naive uniaxial-swelling estimate \((8/\sqrt{\pi})\gamma=1.128\) is identified as incorrect: it misses the outer elastic response to the membrane force in the constrained layer. Loading-2 numerics recover \(\lim\sqrt{s}\cdot s\tilde V_2=-8\eta(1-\nu_u)=-1.6020\) to about \(10^{-4}\) (`ANSWER.md` 315–317; independent check: exact \(-1.602\)).

**Tips (same order).** The undrained tip suction \(p_K=-2\beta K\cos(\theta/2)/\sqrt{2\pi r}\) is harmonic and has \(\zeta=0\), but the face flux is non-integrable. For \(t>0\) it is regularised in \(r\sim\sqrt{ct}\). First-order (mechanics frozen): \(p=p_K\,\mathrm{erf}(r/2\sqrt{ct})\), checked symbolically (`check_erf_solution.py`; PDE residual 0, vanishes on \(\theta=\pm\pi\)). That gives

\[
A_{\mathrm{tip}}^{(1)}=4\sqrt{\pi}\,\gamma(1-\nu_u)=1.06347.
\]

Full coupling: \(A_{\mathrm{tip}}=(2/\sqrt{\pi})(1-\nu_u)\hat{c}(\gamma)\) with \(\hat{c}(0.25)=1.6043\), so \(A_{\mathrm{tip}}=1.0862\). Reported \(A=1.3541+1.0862=2.440\) (and 2.44 in `results.json`) against key 2.440220.

The first-order closed form \(A\approx 4(4+\pi)\gamma(1-\nu_u)/\sqrt{\pi}=2.4175\) is disclosed as 0.9% low. It is a weak-coupling approximation, not the reported answer.

Leading remainder \(O((ct/L^2)^{3/4})\) is consistent with the \(t^{1/4}\) correction to \(K\) feeding the volume at the next order.

---

## 4. Part (b) — correct

| Item | Verdict | Reason |
|---|---|---|
| `K2_0=0` | correct | Inner problem \(O(t^{1/4})\); paper Eq. 38. Reported 0. |
| `K2_inf=-\eta` | correct | Coupled 2-D back-stress, not \(-\alpha\). Reported \(-0.33375\). |
| `V2_inf=-\pi\alpha(1-2\nu)/2\) | correct | Closed form, also \(-\pi\eta(1-\nu)\). Reciprocal check using Griffith integrals \(\int\sigma_{kk}\,dA=-\pi L^2\), \(\int\sigma_{kk}^2\,dA=2\pi L^2\) (`griffith_integrals.py`). Reported \(-0.838805\). |

The logarithmic approach and the small-time \(K_2\sim -0.2599\,\tau^{1/4}\) are supplementary (see §7). They support the limits but are not requested answers.

Linear superposition is allowed to close the crack (\(V_2<0\)); the submission notes the contact caveat (`ANSWER.md` 341–342). That is a modelling remark, not an error.

---

## 5. Part (c) and the nine numbers — correct

Independent evaluation of the closed forms and the 1% / 0.01 tests (`scratch/independent_checks.py`):

| Quantity | Reported | Reference | Error | Verdict |
|---|---:|---:|---:|---|
| `K1_0` | 0.80759 | 0.807593 | rel \(3.7\times 10^{-6}\) | correct |
| `K1_inf` | 1.0 | 1 | 0 | correct |
| `V1_0` | 1.884956 | \(\pi(1-\nu_u)\) | rel \(2\times 10^{-7}\) | correct |
| `V1_inf` | 2.513274 | \(\pi(1-\nu)\) | rel \(5\times 10^{-8}\) | correct |
| `V1_short_coeff` | 2.44 | 2.440220 | rel \(9\times 10^{-5}\) | correct |
| `p_faces_0` | 0.749064 | \(2\beta=0.74906367\) | rel \(4\times 10^{-7}\) | correct |
| `K2_0` | 0.0 | 0 | abs 0 | correct |
| `K2_inf` | \(-0.33375\) | \(-\eta\) | 0 | correct |
| `V2_inf` | \(-0.838805\) | \(-\pi\alpha(1-2\nu)/2\) | rel \(3\times 10^{-7}\) | correct |

Principle checks on the agent’s own numbers: \(V_1(\infty)/V_1(0^+)=(1-\nu)/(1-\nu_u)\); \(K_1(\infty)=1\).

`p_faces_0` is the undrained pressure *in the solid next to the faces*, before diffusion, while the faces themselves have \(p=0\). That is what the problem asks for.

---

## 6. Solver, residuals, convergence, extrapolation

### Formulation

Laplace transform; even finite crack; undrained Papkovich–Neuber field (\(\nu_u\), \(\zeta=0\), \(p\) harmonic) plus one irrotational diffusion mode \(u=\nabla\phi\), \(\nabla^2\phi=(s/c)\phi\) (`ANSWER.md` 231–247). Galerkin basis \(\sqrt{1-x^2}\,U_n\) for opening and \(T_n/\sqrt{1-x^2}\) for flux: the drained \(r^{1/2}\) opening and \(r^{-1/2}\) tip flux. Drained/large-\(k\) integrals in closed form (Weber–Schafheitlin); oscillatory remainder by Gauss–Legendre; Debye tail for \(k>k_{\max}\); \((0,0)\) flux entry by Neumann’s product formula. First solver (`crack_solver.py`) lacked the Debye tail and showed \(10^{-4}\) noise in \(K\); production code is `crack_solver2.py`.

The kernel algebra in `symbols.py` was not re-derived here (that script pickles to a path outside this workspace). Limits of the implemented remainders were checked: \(s\to\infty\) undrained, \(m\to|k|\) drained with \((1-\gamma)/(1-\nu_u)=1/(1-\nu)\). Several global identities would fail if the coupling blocks were wrong (energy, face-layer coefficient, \(\gamma\)-collapse, small-\(s\) loading-1 limits). That is evidence, not a line-by-line proof of the Fourier symbols.

### Saved residuals and convergence

- `large_s_1e+07.json`: \(s=10^7\), \(N=150\) vs \(190\): \(sK_1\) agrees to \(4\times 10^{-9}\); \(|b_{\mathrm{last}}|\) drops from \(4.6\times 10^{-6}\) to \(2.9\times 10^{-7}\).
- `large_s_1e+06.json`: \(N=80,110,140\), \(sK_1\) stable at \(10^{-8}\).
- `conv_test2_large_s.txt`: doubling \(k_{\max}\) changes \(sK_1\) at the \(10^{-8}\) level once \(N\) is large enough that \(|b_{\mathrm{last}}|\) is small. At \(s=10^5\), \(N=20\) is unconverged (\(|b_{\mathrm{last}}|=0.44\)); \(N=60\) is converged. They did not use the unconverged points.
- Small \(s\): graded \(k\)-panels. An ungraded version produced a wrong long-time \(K_2\); that bug was caught because it contradicted the disc argument (Sec. 3). `small_s.log` at \(s=10^{-12}\): \(sK_1=1\), \(sV_1=\pi(1-\nu)\) to the printed digits; \(K_2\) still approaching \(-\eta\) logarithmically, as it must.

### Extrapolation for `K1_0` and `V1_short_coeff`

`final_extrapolation.py` uses the initial-value theorem:

\[
K_1(0^+)=\lim_{s\to\infty}s\tilde K_1,\qquad
A=\frac{2}{\sqrt{\pi}}\lim_{s\to\infty}(s\tilde V_1-V_u)\sqrt{s}.
\]

The clean quantity is \(W=s(\tilde K_1+2\beta\tilde K_2)\), i.e. faces loaded with the undrained pressure \(p=2\beta\), so there is no face layer. Independent re-fit of the saved table (\(s=10^3\ldots 10^7\), powers \(s^0,s^{-1/2},s^{-1}\)):

- \(a=0.807592986\) (claimed 0.807592986);
- max residual \(1.6\times 10^{-8}\);
- \(W(s=10^7)-a=8.2\times 10^{-6}\), which is the last \(s^{-1/2}\) correction, not the error.

Claimed uncertainty \(a=0.807593\pm 5\times 10^{-6}\) is justified as an extrapolation estimate if the \(s^{-1/2}\) model (theoretically expected once the face layer is removed) is accepted. It is not a rigorous bound from a remainder theorem. Rounding the reported `K1_0` to 0.80759 is appropriate.

For \(A\): direct extrapolation 2.43999; exact face + numerical tip 2.44022; numeric face coefficient 1.35382 vs exact 1.354055 (difference \(2.3\times 10^{-4}\)). The slower \(s^{-1/4}\) contamination in \(V\) is real. Claimed \(A=2.4401\pm 0.0003\) is consistent with that spread. Reporting `V1_short_coeff=2.44` matches the uncertainty.

### Representative reruns (copied solver, `scratch/`)

The full \(s=10^7\), \(N=190\) campaign and the Talbot table were **not** rerun.

| Check | Result |
|---|---|
| \(s=100\), \(N=20\) | bit-identical to `large_s_1e+02.json` |
| \(s=1\), \(N=24\), energy identity \(J_{\mathrm{tip}}=V/s+\tfrac12 dV/ds\) | relative residual \(4\times 10^{-10}\) |
| \(s=10^{-4}\), \(N=16\) | matches `small_s.log` line 1 |

Talbot histories at \(\tau=10^{-4}\): \(K_1=0.8273\) vs inner-law \(0.8076+0.1947\tau^{1/4}=0.8271\); \(K_2=-0.0260\) vs \(-0.2599\tau^{1/4}\). At \(\tau=10^4\): \(K_2=-0.3049\) vs log law \(-0.3045\). Short-time \(V_1\) at \(\tau=10^{-4}\) is 1.9065 vs \(V_1(0^+)+2.44\sqrt{\tau}=1.9094\); the \(O(\tau^{3/4})\) remainder is visible and of the expected sign and size. That is a consistency check, not an independent determination of \(A\).

---

## 7. Formulas vs fitted or conjectural extras

**Requested formulas that are closed-form and correct**

- `K1_inf` \(=\sigma\sqrt{\pi L}\)
- `V1_0` \(=\pi(1-\nu_u)\sigma L^2/G\)
- `V1_inf` \(=\pi(1-\nu)\sigma L^2/G\)
- `p_faces_0` \(=2(\nu_u-\nu)\sigma/[\alpha(1-2\nu)]\)
- `K2_0` \(=0\)
- `K2_inf` \(=-\eta\sigma\sqrt{\pi L}\)
- `V2_inf` \(=-\pi\alpha(1-2\nu)p L^2/(2G)\)
- `V1_short_coeff` face piece \(\frac{16}{\sqrt{\pi}}\gamma(1-\nu_u)\) exact
- definitions of \(\eta,\beta,\gamma\)

**Requested items that are numerical, with a justified general definition**

- `K1_0` \(=a(\gamma)\sigma\sqrt{\pi L}\), \(a\) the drained-tip / undrained-far-field factor. Weak-coupling \(a=1-3\gamma/4+O(\gamma^2)\) is derived. No elementary closed form is claimed.
- `A_tip` through \(\hat{c}(\gamma)\); first-order \(\hat{c}=2\pi\gamma\) is derived and matches \(\gamma\to 0\).

**Supplementary / fitted / not requested** (not used as the graded answers)

| Claim | Status |
|---|---|
| Polynomial \(a\approx 1-0.75\gamma-0.0634\gamma^2-0.0582\gamma^3\) | Fit to tabulated \(a(\gamma)\). Max \(\lvert\Delta a\rvert=8.1\times 10^{-5}\) on \(\gamma\in(0,0.5]\), inside the claimed \(10^{-4}\). Not a theory. |
| Polynomial \(\hat{c}\approx 2\pi\gamma+0.434\gamma^2+0.37\gamma^3\) | Fit. Max \(\lvert\Delta\hat{c}\rvert=7.7\times 10^{-4}\) (at \(\gamma=0.5\)); claimed \(\lvert\mathrm{error}\rvert<5\times 10^{-4}\) is **slightly overstated**. Does not affect \(A\) at \(\gamma=0.25\) beyond the already-quoted \(\pm 3\times 10^{-4}\). |
| \(K_1\approx 0.8076+0.1947\tau^{1/4}\) | Follows from \(K_1(t)=a K_\infty-2\beta K_2(t)\) plus the loading-2 inner scaling; coefficient from the same large-\(s\) fit. Supplementary. |
| \(K_2\sim -0.2599\,\tau^{1/4}\), Laplace \(s^{1/4}s\tilde K_2\to-0.23554\) | Numerical inner-layer coefficient. Consistent with \(\Gamma(5/4)\). Supplementary. |
| Log laws with \(D=4\ln 2-2\gamma_E\) | The constant matches the 2-D logarithmic capacity of a segment; \(C=-X_\infty\) from a three-point fit at tiny \(s\). Strongly motivated, still a fit of the prefactor. Supplementary. |
| Bounds \((1-\gamma)^{3/4}<a<\sqrt{(1-\gamma)(1-\gamma/2)}\) | Hold on the tabulated points; the upper bound is the weak-coupling \(J\) result. Not a proof for all \(\gamma\). |
| Yu-type \([(1-\nu_u)/(1-\nu)]^{0.743}\) | **Not claimed.** At \(\gamma=0.25\) it gives 0.80755, within \(5\times 10^{-5}\) of \(a\); at \(\gamma=0.5\) it is 0.7% low. The submission was right not to treat it as a general formula. |
| Combined fluid-filled crack \(K\to(1-\eta)p_f\sqrt{\pi L}\) | Immediate superposition; not requested. |

General formulas in `results.json` are therefore acceptable as the requested strings: the closed forms are exact, and `K1_0` / `A_tip` are stated as inner-problem definitions plus numerics, not as fake closed forms.

---

## 8. Verification coverage

Covered by the submission, with independent confirmation of a subset:

1. Analytic limits \(s\to 0\) (loading 1) and \(s\to\infty\) (\(s\tilde V_1\to\pi(1-\nu_u)\)).
2. Laplace-domain energy / \(J\) identity, including the fluid tip term (rerun at \(s=1\)).
3. Mode and quadrature convergence, including a documented failure mode (missing Debye tail; ungraded small-\(s\) panels).
4. Single-parameter reduction to \(\gamma\).
5. Weak-coupling slopes as \(\gamma\to 0\).
6. Face-layer coefficient for loading 2 against the exact 1-D result.
7. Long-time 2-D log fit, capacity constant.
8. Griffith integrals used in \(V_2(\infty)\).
9. Symbolic erf tip field.
10. Talbot histories vs small-time and log laws.

Not covered, and not assumed:

- Closed Wiener–Hopf factorization of \(a(\gamma)\).
- Independent re-derivation of every Fourier-symbol algebra line.
- Rerun of \(s=10^6\)–\(10^7\) at \(N\sim 190\), or the \(M=12\) Talbot campaign.
- A 2-D time-domain finite-element check of the same BVP (the key cites one; it was not repeated here).
- Contact / finite deformation / nonlinear permeability.

Numerical agreement with the nine key values is a **consistency** test. The case that the discrete system is the right BVP rests on the face conditions, the drained-tip extraction, the identities, and the representative reruns, not on the numbers alone.

---

## 9. Reference and validation limits

- The grading key is `REFERENCE.md` as corrected in log P374, not the raw Detournay–Cheng 1991 numbers for `K1_0` and `V1_short_coeff`.
- Paper values 0.75 and 1.354 are recorded so that matching the paper against the corrected key can be classified as wrong. This submission does not match those two paper values.
- The problem statement’s \(t=0^+\) remains order-of-limits ambiguous. This grade does **not** treat an undrained-outer answer of 1 as incorrect; it treats it as partial if it is the only `K1_0` given. Submission C gives the drained-tip value and documents the other reading.
- Tolerances: 1% relative (nonzero), 0.01 absolute (`K2_0`). All nine pass by a wide margin.
- Independent computation in this review: closed-form arithmetic; re-extrapolation of *saved* large-\(s\) data; polynomial-fit residuals; three solver reruns at modest \((s,N)\); energy identity at \(s=1\). Not a reproduction of the whole campaign.
- `write_results.py` and `symbols.py` contain absolute paths to another workspace. That is a provenance artifact. The `results.json` sitting in `C/` is consistent with `final_extrapolation.json` and with the closed forms.

---

## 10. Verdicts by part

| Part | Verdict | Summary |
|---|---|---|
| (a) | **correct** | Undrained/drained volumes and \(K_1(\infty)\) exact. `K1_0` is the drained-tip inner limit, with the outer reading recorded. Short-time \(A\) includes faces (exact) and tips (numerical), not the paper’s faces-only coefficient. |
| (b) | **correct** | \(K_2(0^+)=0\); coupled 2-D back-stress \(-\eta\) and \(V_2(\infty)\) with the nonuniform logarithmic far-field structure. Sign of loading 2 correct. |
| (c) | **correct** | All nine values inside tolerance; formulas match the derivations; extras are labelled. |
| Verification | **sufficient** for the numerical claims, not a closed-form existence proof | Residuals, convergence, identities, and representative reruns support the BVP. Uncertainty on \(a\) is a well-behaved extrapolation estimate; on \(A\) it is limited by \(s^{-1/4}\). The \(\hat{c}(\gamma)\) fit error bound is slightly tight. |
| **Overall** | **correct** | |

No requested item is missing. No core answer is unsupported. Fitted polynomials, small-time \(t^{1/4}\) prefactors, and log-capacity constants are supplementary and were graded as such.
