# Key for problem B: a suddenly pressurized crack in a poroelastic solid

**Source:** E. Detournay & A. H.-D. Cheng 1991, "Plane strain analysis of a stationary hydraulic
fracture in a poroelastic medium", *Int. J. Solids Struct.* 27(13), 1645–1662, doi:10.1016/0020-7683(91)90067-P.

**Status:** written 8 October 2026 (log P341), corrected on 25 and 26 October 2026; `../FIXES.md` records the fixes and
Teng's decisions of 6 November 2026. Open book; every element is cited. Our own
computations are marked "ours".

## Key correction (log P374): K1_0 and V1_short_coeff

For these two quantities the paper gives a simplified solution (K1_0) and the face term alone (V1_short_coeff).
The first key took those values as exact; the key now uses the full solution, and the paper's values are kept on
record below.

| Quantity | Paper's simplified solution | Key (full solution) | Evidence |
|---|---:|---:|---|
| K1_0 | 0.75 (Eq. 29, "heuristic argument") | **0.807593** | GPT MP-37/38/39; Codex (Wilson factorization, 0.8075930); Grok (0.80759 ± 0.00005); our FE (0.80 ± 0.01) |
| V1_short_coeff | 1.354055 (Eq. 18, faces only) | **2.440220** | 1.354055 (faces) + 1.086165 (tips): Codex and Grok; GPT MP-37/38; our FE (2.2–2.4, extrapolated) |

- **Literature:** for gels (ν_u = ½), Yu et al. 2020, Eq. 12 (from Atkinson & Craster) write the tip limit as
  K(0⁺) = K_∞[(1 − ν_u)/(1 − ν)]^α with α ≈ 0.735 (their α, not the Biot coefficient). For this material the
  factorization corresponds to an exponent 0.743 (ours). Details: `CASE_STUDY_GPT_CHALLENGE.md`, `FE_CHECK.md`.
- **Order of limits.** The problem says only "just after loading (t = 0⁺)". GPT showed this is
  ambiguous:
  - the tip limit (r → 0 at fixed t > 0, then t → 0⁺) gives 0.807593, the key;
  - the reverse order gives exactly 1, the undrained outer intensity (Codex);
  - 0.75, the paper's simplified solution, is neither limit.
- **Grading:**
  - an answer of 1 is "partly": right under one reading of the statement, but without the drained
    tip, which is the decisive point;
  - an answer of 0.75 is not the value asked for: it is the paper's simplified solution; record that;
  - V1_short_coeff = 1.354 is not the value asked for: it is the face term alone, without the tips.
- **For the next round,** the statement fixes the order of limits (as in G1 and G2).

## Sources for each element of the problem

| Element | Paper |
|---|---|
| Infinite, homogeneous, isotropic Biot medium; plane strain; stationary crack of length 2L | p. 1645 (assumptions); Fig. 1 |
| Constitutive Eqs., equilibrium, Darcy, mass balance; c | Eqs. 1–5, 8 (pp. 1646–1647) |
| The two loadings and zero initial conditions | Eqs. 9–10 (p. 1648) |
| Crack volume and K as the reported quantities | pp. 1648–1652 |
| ν = 0.2, ν_u = 0.4, α = 0.89 | the paper's examples (p. 1655: "The material constants used for these examples are ν = 0.2, ν_u = 0.4 and α = 0.89") |

## Reference answer (the paper's asymptotic results, pp. 1648–1652, with the full solution for K(0⁺) and the small-time volume)

Here η = α(1 − 2ν)/[2(1 − ν)] (Eq. A3, p. 1661; for (c), η = 0.33375, ours).

**(a) Loading 1.**

*At t = 0⁺ the solid responds with its undrained moduli.*
- **Opening:** the Sneddon ellipse with ν_u, D_n = 2(1 − ν_u)√(L² − x²)/G (Eq. 13, p. 1648).
- **Volume:** V(0⁺) = πL²(1 − ν_u)/G (Eq. 15, p. 1649).
- **Induced pore pressure** on the faces: p = 2(ν_u − ν)/[α(1 − 2ν)] (Eq. 17, p. 1649).

*Stress intensity at t = 0⁺:* **K(0⁺) = 0.807593 √(πL)** (tip limit r → 0 first, then t → 0⁺), not √(πL).
- **Why the tip is drained (p. 1650):** the induced pore pressure has "infinite gradients" at the tips,
  so "the crack tips are effectively drained at t > 0". The r^−½ singularity is preserved (Rice &
  Simons 1976; Simons 1977).
- **Full solution:** the drained tip field is matched to the undrained crack-scale field through the
  diffusion layer, a 2 × 2 Wiener–Hopf problem whose symbol depends only on η_t = (ν_u − ν)/(1 − ν_u),
  here 1/3. Its factorization gives 0.807593 (key correction above).
- **The paper's simplified solution:** it puts the undrained opening into the drained opening–K relation,
  K = G√π / [√8 (1 − ν)] · lim D_n/√ξ (Eq. 27), and gets K(0⁺) = [(1 − ν_u)/(1 − ν)] √(πL) = 0.75√(πL)
  (Eq. 29, p. 1651), which leaves out the matching layer. The paper calls it a "heuristic argument"; its
  numerical solution approaches that value at small times (Fig. 9, p. 1659).
- **This is the decisive point of the problem:** a naive "undrained elastic" answer gives √(πL).

*Long times:*
- V(∞) = πL²(1 − ν)/G (Eq. 19, p. 1649), so the volume grows by πL²(ν_u − ν)/G (Eq. 20).
- K(∞) = √(πL) (Eq. 26, p. 1650), "independent of any material property".

*Small times:* V(t) = V(0⁺) + A · (L²/G) √(ct/L²) with **A = 2.440220** for (c)'s values (key correction above).
- **Faces (the paper's Eq. 18, p. 1649):** 16(1 − ν_u)(ν_u − ν) / [√π(1 − ν)] = 1.354055. The route:
  1. the induced pore pressure (Eq. 17) is dissipated by one-dimensional diffusion through the faces
     (Eq. 21);
  2. the resulting leak-off (Eq. 22) is equal, by the reciprocal theorem (Eqs. 30–31), to the
     loading-2 volume change;
  3. that change, scaled by the induced pressure, gives Eq. 18.
- **Tips:** near each tip the drained region shrinks like √(ct), but the stress trace of the undrained slit
  solution is singular there, so the tips add a term of the same √t order: 1.086165, from the same
  factorization (ours; not checked independently, our FE extrapolates A to 2.2–2.4).
  *Note, 6 November 2026:* a direct finite-crack solve on 29 October 2026 gave the tip term as 1.09 ± 0.01
  (`../independent_check_2026-10-29/REPORT.md`); the key keeps 1.086165 (`../FIXES.md`, decisions of 6 November).
- The paper calls its small-time result "an excellent approximation for t* < 0.001" (p. 1656, Fig. 6).

**(b) Loading 2.**
- K(0⁺) = 0 (Eq. 38, p. 1652): "there is no displacement on the crack face induced by mode 2 loading".
- K(∞) = −η√(πL) (Eq. 39, p. 1652). At long times a uniform confining stress σ_xx = σ_yy = −η
  develops (Eq. 35), equivalent to a "back-stress" η on the faces.
- V(∞) = −απL²(1 − 2ν)/(2G) (Eq. 33, p. 1651): a closure.

**(c) Numbers (ours, from the paper's formulas; checked against its Figs. 6, 7, 9, 10).**

| Quantity | Value | Figure check |
|---|---:|---|
| K1_0 | 0.807593 (P374; the paper's simplified solution gives 0.75) | Fig. 9 dashed line ≈ 1.33 = 0.75√π (the simplified solution) |
| K1_inf | 1 | Fig. 9 ≈ 1.77 = √π |
| V1_0 | 1.88496 (= 0.6π) | Fig. 6 starts at ≈ 1.89 |
| V1_inf | 2.51327 (= 0.8π) | Fig. 6 dashed ≈ 2.51 |
| V1_short_coeff | 2.44022 (P374; faces 1.35406, the paper's Eq. 18, + tips 1.08617) | Fig. 6 at t* = 0.01: 1.885 + 0.135 = 2.02 (faces only) |
| p_faces_0 | 0.74906 | — |
| K2_0 | 0 | Fig. 10 starts at 0 |
| K2_inf | −0.33375 | Fig. 10 dashed ≈ 0.59 = 0.33375√π |
| V2_inf | −0.83881 | Fig. 7 dashed ≈ 0.84 |

## Grading (our choices)

**Automatic layer** (`checks.py`):
- the nine values: within 1% where the value is nonzero, except K1_0 within 2% and V1_short_coeff within 5%
  (26 Oct: a resolved FE route gives 0.798–0.801 and 2.2–2.4); K2_0 within 0.01;
- **principle checks on the agent's own numbers** (independent of the key):
  - V1_inf/V1_0 = (1 − ν)/(1 − ν_u). Both volumes are Sneddon ellipses with the drained and undrained
    ν.
  - K1_inf = 1 in units of √(πL), the Griffith value: the drained long-time limit carries no material
    constant.

**Grader layer:**
- whether K1_0 is reasoned with drained tips or simply taken as the elastic value;
- the reasoning for the small-time growth (diffusion, reciprocity);
- the sign of loading 2;
- every error, classified; checks able to fail; source naming.

## Lessons this problem probes

- **L10 (knowing the boundary):** the drained crack-tip zone at t = 0⁺.
- **L9 (conventions):** ν against ν_u; the definition of V; signs.
- **L1:** closed forms exist.
- **L3:** the small-time growth needs the right method (1D diffusion plus reciprocity), not a full
  transient solution.
- **L8:** source naming.
