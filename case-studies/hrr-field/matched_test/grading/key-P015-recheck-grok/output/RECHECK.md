# Recheck: crack_tip_asymptotics v0.2

Started: 02:22:29 UTC. Complete: 02:46:42 UTC.

Scripts: `scratch/recheck_probes.py`, `scratch/more_probes.py`, `scratch/existence_table.py`,
`scratch/modeII_double_root.py`, `scratch/crash_at_table_angle.py` (copies of the library in
`scratch/crack_tip_asymptotics/`; `./library/` not edited).

## Checklist

- [x] 1. Previous-review findings: resolved / partly / not resolved / rejected
- [x] 2. New TOPIC §6 and §3 existence table
- [x] 3. New problems from the revision (break new checks)
- [x] Verdict

## 1. Previous findings

Reran the v0.1 breakers on the revised `crack_tip_checks.py`.

### Mechanics / topic (PREVIOUS_REVIEW §1)

| Finding | Status | Evidence |
|---|---|---|
| Mode II Williams equation missing | **resolved** | TOPIC.md:67; `williams_exponent(..., mode="II")`; 150° notch p = 0.269099 (`scratch/recheck_probes.py`, matches the old independent value). |
| Polar shear 2ε_rθ vs tensor ε_rθ | **resolved** | TOPIC.md:50 and checker docstring state tensor shear; displacement identity still `Ur' + q Ut − 2 ert` (`crack_tip_checks.py:291`). |
| SSY K replaces α by α/(1−ν²) | **resolved** (as documentation) | ENTRY.md:31, 41: plane-strain (1−ν²) is in the remote elastic J only. No new claim. |
| Printed vs computed Eₙ ≤ 0.85% | **resolved** (still flagged) | ENTRY.md:30 and `reference_data.json` unchanged; still do not treat printed Eₙ as 12-digit truth. |
| Airy signs, incompressibility, Cartesian ∂u/∂x₁ | **rejected** as issues | Unchanged and still consistent (TOPIC §2, §4). |

### False positives (correct fields fail/crash)

| Finding | Status | Evidence |
|---|---|---|
| Relative Cartesian / Saint-Venant residuals reject T-stress and mode II | **resolved** | Residuals / field (or strain) scale at r = 1 (`crack_tip_checks.py:169–171, 218, 234`). T-stress `leading=False`: Cartesian 2.35e-10 PASS, Saint-Venant 3.10e-7 PASS. Mode II crack: Cartesian 9.11e-9, Saint-Venant 3.05e-7 PASS (`scratch/recheck_probes.py`). |
| 90° wedge: `williams_exponent` raises; `check_field` crashes | **resolved** for mode I | `williams_exponent(π/2, "I")` = 0; uniform σ_yy field PASSes (`test_90_degree_wedge…`; probes). Field is the Williams p = 0 Airy F = 1 + cos 2θ (`scratch/more_probes.py`). |
| Mode I only for J / displacement; p = 1 divides by 1+q and crashes | **resolved** (no crash) | `1+q ≈ 0` returns n/a (`crack_tip_checks.py:273–277`). Williams stresses at p = 1: no exception; FAIL energy + eq + crack exponent (`scratch/recheck_probes.py`). J/displacement still mode I only (documented TOPIC.md:151–152). |
| Coarse 21-point full-range Williams FAIL | **partly** | Still FAILs angular eq (5.7e-3) and Cartesian (2.2e-3). Now emits a resolution n/a (`crack_tip_checks.py:179–182`). Half-range 21 pts also FAILs (mirrors to 41). Documented “refine first”; still a trap. |

### False negatives / leading-term

| Finding | Status | Evidence |
|---|---|---|
| Higher-order Williams (λ = 1.5, p = −0.5) PASSes except crack exponent | **resolved** via `leading` | `leading=True`: FAIL crack exponent (and a non-singular report). `leading=False`: all PASS (`scratch/recheck_probes.py`). |
| Linear notch Δp = 10⁻⁵ FAIL Williams, eq PASS; power-law notch no p-oracle | **partly** | Linear 120°: Δp = 10⁻⁴ FAILs angular eq (1e-4 > tol_eq 1e-5) **and** Williams; Δp = 10⁻⁵ FAILs Williams only, angular eq = 1e-5 PASS (on the tolerance). Power-law Δp = 10⁻⁴ is now intended to fail angular eq (`test_power_law_notch_exponent_error_of_1e4_is_caught`; not re-run here: needs frozen HRR paths). |
| Negative-K_I closing crack PASSes | **resolved** | Opening check: U_θ(π) = +1.452 FAIL only that check; +K_I PASSes (`scratch/recheck_probes.py`). |
| Small T mixed into Williams at p = 1/2 | **partly** | ε = 10⁻⁴ now FAILs **angular** eq (1.15e-4 > 1e-5); Cartesian still PASS. Tighter than v0.1 (then ε ≤ 10⁻⁴ passed). Residual pollution of a leading term is still not a dedicated check. |
| F = cos³(θ/2) is Williams, not a mutant | **rejected** as a bug | Unchanged. Fake F = cos⁵(θ/2) and (1+cos θ)² still FAIL compatibility only (`scratch/recheck_probes.py`). |

### Solver constraints (PREVIOUS_REVIEW §3)

| Finding | Status | Evidence |
|---|---|---|
| J/displacement advertised but mode I only | **resolved** (documented) | TOPIC.md:151–152; mode II returns n/a (`crack_tip_checks.py:267–270`). Mode II crack still has no J number. |
| Check 7 forces the leading J exponent | **resolved** via `leading=False` | Documented TOPIC.md:148–149 and checker docstring. Default `leading=True` on T-stress still FAILs “crack exponent fixed by J” **and** emits the non-singular report (`scratch/recheck_probes.py`). Easy to misuse if the flag is left at default. |
| `test_checks.py` hard-coded v1 / T2 paths | **not resolved** | `test_checks.py:17–21` still imports `PILOTS/corpus/knowledge_v1/...` and `<probe-runs>/P015-build-codex/T2/reference.py`. |
| hrr.py crack-only; α name clash | **resolved** (documented) | ENTRY.md:82–86; TOPIC now uses β for the wedge (`TOPIC.md:26–27`). |
| Dimensional constants = 1 | **resolved** (stated) | TOPIC.md:100; checker docstring. |

### Missing for related problems (PREVIOUS_REVIEW §4)

Mode II characteristic equation: **resolved**. 90° / p = 0: **resolved** for mode I. Relative residuals for T-stress: **resolved**. Contact/opening: **resolved** for mode I cracks. Interface / flow / 3D: still out of scope (TOPIC §5–6); that is appropriate.

Minimum v0.1 fixes (absolute residual floor; p = 0 Williams; document leading/J/mode I): **done**.

## 2. New section 6 and the §3 existence table

### Existence table (TOPIC.md:71–78)

Computed with `williams_exponent` (`scratch/existence_table.py`):

| 2β | table p_I | computed p_I | table p_II | computed p_II |
|---|---|---|---|---|
| 180° | 0 | 0 | none singular | −1.0 (λ = 2; bounded) |
| 240° | 0.384 | 0.384269 | none singular | −0.1489 (bounded) |
| 257.45° | 0.431 | 0.431701 | 0 (double root) | −2.6×10⁻⁵ at 257.45°; **crash at the exact root** |
| 270° | 0.456 | 0.455516 | 0.092 | 0.091471 |
| 300° | 0.488 | 0.487779 | 0.269 | 0.269099 |
| 360° | 0.5 | 0.5 | 0.5 | 0.5 |

Physics of the rows is right: mode I singular only for a re-entrant corner (p_I(179°) < 0 < p_I(181°)); mode II only past the double root of tan 2β = 2β at 2β = 257.45340° (`scratch/existence_table.py`, `test_singularity_thresholds`). “None singular” at 180°/240° mode II is honest (p < 0). Labels: **[derived]** is fair for the numbers (they come from this topic’s root finder). Rounding 0.4555→0.456, 0.09147→0.092 is fine.

**Overstated / hole:** the table cell “0 (double root)” at 257.45° is the physics, but `williams_exponent(β_crit, "II")` **raises** `ValueError: no real root` because λ = 1 is stripped for all mode II (`crack_tip_checks.py:66–69`), including the isolated angle where that root **is** the stress-free double root. Crash window ~0.005° above the critical angle (`scratch/modeII_double_root.py`). `check_field` at 257.45° (slightly below) does not crash (`scratch/crash_at_table_angle.py`). The exception named in TOPIC.md:67–68 is not implemented.

### Section 6 rows

| # | Physics | Label honesty | Notes |
|---|---|---|---|
| 1 Inner limit | Right: small-strain singularity is inconsistent at r → 0; process zone / blunting. | “a few δ_t”, δ_t ≈ d_n J/σ₀, d_n ≈ 0.3–0.8 **[from memory: McMeeking 1977; Shih 1981]** — correctly labelled. | Nothing overstated. |
| 2 Outer limit | Right: T / Q / A₂; J-dominance needs ligament ≫ J/σ₀. | Larsson & Carlsson **[standard]** OK; size factors 25 J/σ₀ bending / 200 J/σ₀ tension **[from memory]** OK. | Fine. |
| 3 Singularity exists | Right, matches §3. | **[derived]** OK. | “Failure is a strength question” is the right alternative. |
| 4 Distinct eigenvalues | Right: double root → r^λ ln r; leading term ill-conditioned nearby. | Double root **[derived]**; Sternberg–Koiter at the same angle **[from memory]** — honest. | The link is standard lore; not checked against Sternberg–Koiter here. |
| 5 Faces open | Right: traction-free assumes separation; U_θ(π) > 0 is closing. | Opening **[derived, check]**; Comninou contact **[from memory]** OK. | Sign convention matches the check (U_θ(β) < 0 opens). Interface interpenetration is from memory, appropriate. |
| 6 Proportional loading | Right: inside the separable zone σ ∝ J^{1/(n+1)} × fixed tensor. | That piece **[derived]**; unloading/growth/cyclic **[standard]** OK. | ENTRY.md:70–72 still warns Hutchinson did not establish proportionality. |
| 7 Superposition | Right for power-law mixed mode. | Mixity **[from memory: Shih 1974]** OK. | Fine. |
| 8 2D | Right: plane strain near the front; plane-stress HRR needs a zone thicker than the plate. | **[standard]** OK. | Fine. |
| 9 Field vs fracture | Right: K/J as energy parameters vs notch amplitude (stress × length^p). | Leguillon **[from memory, P040 not read]** — explicitly honest. | Fine. |

**Missing (not fatal):** no row on **finite strain / blunting as a separate constitutive change** (item 1 covers geometry of the inner zone; Knowles–Sternberg is in §5). No quantitative annulus for HRR vs T/Q. Mode III / antiplane is absent (out of Hutchinson’s paper). Nothing important is overstated given the labels.

## 3. New problems introduced by the revision

Tried to break: crack opening, non-singular report, mode II Williams, `leading`, field-scale normalization.

1. **Mode II Williams at the table’s critical angle still crashes** — `williams_exponent` and therefore a linear-notch `check_field(..., mode="II")` at the exact double root. New relative to the 90° mode-I fix. Evidence: `scratch/modeII_double_root.py`, `scratch/existence_table.py`. Small: isolated angle; tests use crit ± 0.005 rad and pass.

2. **`leading=True` (default) + p ≤ 0 on a crack:** emits the non-singular report **and** FAILs “crack exponent fixed by J” (T-stress, higher-order p = −0.5). TOPIC §7 says `leading=False` skips those checks; a user who only reads “p ≤ 0 → non-singular report” still gets a FAIL. Evidence: `scratch/recheck_probes.py`. Small documentation / control-flow inconsistency, not a wrong residual.

3. **Opening check is a strict `U_θ(β) < 0` with no tolerance.** T-stress on a crack PASSes opening because U_θ(π) ≈ −7.9×10⁻¹⁶; hydrostatic U_θ = 0 FAILs opening (and faces). A roundoff of the opposite sign would FAIL a valid T-stress. Evidence: `scratch/recheck_probes.py`. Low.

4. **Opening is not applied to notches** (report only) and not to mode II. Closing a 120° notch PASSes everything (`scratch/recheck_probes.py`). Intended in the note; still a silent pass if faces would close.

5. **Zero field / `scale = 0`:** divide-by-zero → NaN residuals → all FAIL, plus RuntimeWarning (`scratch/more_probes.py`). Pathological.

6. **Field-scale residual is sound** for T-stress and mode II (the v0.1 high FPs). Tiny (1e-16) Williams still PASSes. Mixing a regular term into a leading p = 1/2 field is now easier to fail on angular eq because tol_eq = 10⁻⁵ — a feature for exponent errors, a bit sharp for two-term pollution.

7. **`reported_J` on a notch** still never FAILs (1e9 → n/a). Documented; unchanged.

8. **Normalization report** (max Se of the given angular field) works; Williams K = 1 gives max Se ≈ 0.364, not 1 — report-only, as intended.

9. **90° mode II** returns p = −1 (λ = 2), not a crash; a constructed Airy field PASSes with a non-singular report (`scratch/crash_at_table_angle.py`). Table “none singular” is right.

No mechanics regression found in Airy identities, J = 1−ν² for K_I = 1 plane strain, or the opening sign vs textbook u_y (`test_crack_face_displacement_matches_textbook_value` in the library tests; not executed here because those tests import frozen v1 paths).

## Verdict

**Freeze after the listed small fixes.**

1. In `williams_exponent` mode II, do **not** drop λ = 1 when it is a double root (tan 2β = 2β); return p = 0, matching the table. Optionally treat a small window around that angle the same way, or stop raising.
2. If `leading=True` and p ≤ 0, skip “crack exponent fixed by J” (the non-singular report already covers it), **or** state in one sentence that T-stress / higher-order terms must pass `leading=False`.
3. Optional: opening check `U_θ(β) < −tol` (or ≤ 0) so T-stress is not sign-of-roundoff.

The v0.1 high-severity FPs (relative residuals, 90° crash, p = 1 crash) are fixed. §6 is physically right and honestly labelled. Do not block on (3) or on the still-hard-coded `test_checks.py` paths.
