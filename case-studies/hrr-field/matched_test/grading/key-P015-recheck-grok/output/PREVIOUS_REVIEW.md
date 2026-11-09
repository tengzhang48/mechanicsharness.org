# Review: crack_tip_asymptotics

Started: 00:01:48 UTC. Complete: 00:26:48 UTC.

## Checklist

- [x] 1. Mechanics correctness of TOPIC.md and ENTRY.md
- [x] 2. Try to break the checks (false negatives / false positives)
- [x] 3. Solver-constraining or misleading content
- [x] 4. What is missing for a related problem
- [x] Verdict

## 1. Mechanics correctness

TOPIC equilibrium, Airy reconstruction, J scaling, energy, and the mode-I Williams characteristic equation are right. ENTRY matches Hutchinson 1968 on the results that matter (p, Iₙ, Eₙ, SSY K, scope). Evidence: polar equilibrium on σ = r^(−p) Σ; Airy φ = A r^(2−p) F; J ∝ r^(1−p(n+1)); paper eqs. (23)–(27), (32), (34), Table 1 p. 21.

**Low — TOPIC.md:59 / crack_tip_checks.py:42.** Williams is given only for mode I: sin(2λα)+λ sin(2α)=0. Mode II is sin(2λα)−λ sin(2α)=0. Not wrong for the stated mode-I check; incomplete if a solver brings a linear mode-II notch. 120° solid half-angle has no singular mode-II root in (0,1]; 150° has p≈0.269 (scratch).

**Low — TOPIC.md:43 vs crack_tip_checks.py:84.** Polar shear is written as engineering strain 2ε_rθ. The checker stores tensor ε_rθ and uses 2ε_rt in the displacement identity (`shear = Ur' + q Ut − 2 ert`, line 237). Internally consistent; a reader mixing the two conventions could mis-code compatibility.

**Low — ENTRY.md:31 / hrr.py:231–238.** SSY K for plane strain replaces α by α/(1−ν²). That matches Hutchinson (32) and the remote elastic J = π(1−ν²)σ∞² in the paper’s non-dimensionalisation. Correct if the same non-dimensional J is used; easy to mis-copy into a dimensional code that already put (1−ν²) in J = (1−ν²)K_I²/E.

**Low — ENTRY.md:30 / reference_data.json:91–92.** Printed vs computed Eₙ still differs by up to 0.85% with no established cause. Iₙ is within 0.27%. The JSON already flags this; do not treat printed Eₙ as 12-digit truth.

No mechanics error found in the Airy signs (TOPIC Σ_rθ = −(1−p)F′ vs hrr srt = (1−s)f′ with s=2−p), plane-strain incompressibility, J = ∫(W n₁ − t·∂u/∂x₁) ds, or the requirement that ∂u/∂x₁ be Cartesian (the I≈5.65 vs 3.86 note is right).

## 2. Breaking the checks

Scripts: `scratch/break_checks.py`, ad-hoc probes in `scratch/` using a copy of `crack_tip_checks.py`. Williams mode I (KI=1): all applicable checks PASS; reported J = 1−ν² (PE) or 1 (PS). Fake Airy F=cos⁵(θ/2), (π²−θ²)², (1+cos θ)²: equilibrium holds, Saint-Venant fails (not FN). Wrong p, scaled shear, wrong n, plane mismatch, reported I=5.65, truncated faces, broken symmetry, p(n+1)=2: caught as in `mutant_detection.py`.

### False positives (correct fields fail or crash)

**High — crack_tip_checks.py:181–206.** Cartesian equilibrium and Saint-Venant use a *relative* residual with `mag = max(|terms|, 1e-300)`. When the true divergence / incompatibility is ~0, roundoff ~10⁻¹⁰ gives relative residual 1.

- Valid T-stress σ_xx=T, p=0 (Airy-equilibrated, traction-free, mode-I even): angular eq 5.6e-9 PASS, Cartesian 1.0 FAIL, Saint-Venant 1.0 FAIL. Absolute FD residuals are ~10⁻¹⁰. Also fails “crack exponent fixed by J” (p=0≠1/2), which is intended only if this is claimed as the *leading* crack term.
- Valid mode-II crack Airy F=sin(3θ/2)+sin(θ/2), p=1/2: faces, symmetry, both equilibria PASS; Saint-Venant 1.0 FAIL. Worst sample is θ=0, where second-derivative mag ~10⁻¹⁰ and rel=1. At 45°, rel~3e-7.
- Uniform hydrostatic p=0: Cartesian rel=1 at some samples from 10⁻¹² roundoff (also not traction-free, so other checks fail for other reasons).

**High — crack_tip_checks.py:39–54, 225–228.** Linear 90° notch (α=π/2): λ=1, p=0 is a valid non-singular eigenfield (faces traction-free). `williams_exponent_symmetric` raises `ValueError: no root in (0, 1]` instead of returning p=0. `check_field` crashes on a correct field.

**Medium — crack_tip_checks.py:103–109, 230–236.** Displacement/J path is mode I only (`Ut(0)=0`). Mode II never gets a J or displacement-shear check. `p=1` linear (q=−1) divides by `1+q` and can raise `ValueError: y must contain only finite values` (observed at p=1.0 on Williams stresses). Energy already fails at p=1, but a crash is not a clean FAIL.

**Low — coarse grid.** Full-range Williams crack, 21 θ-points: angular eq 5.7e-3, Cartesian 2.0e-2, Saint-Venant 3.9e-2 FAIL. Half-range 21 points (mirrored to 41) PASS. 120° notch, 21 points on [0,α]: PASS. Matches the documented “refine before blaming the model”; still a trap for a solver that samples one side densely and the full interval coarsely.

### False negatives (wrong fields pass)

**Medium — leading-term vs any eigenfield.** Higher-order Williams crack mode (λ=1.5, p=−0.5): faces, eq, compatibility, energy all PASS; only “crack exponent fixed by J” FAILs. Intended if the claim is the J-leading term; a wrong *subdominant* field with a faked p=1/2 would still have to satisfy equilibrium at that p (angular residual ~O(Δp)), so this is not a free pass for a wrong leading field.

**Medium — crack_tip_checks.py:216–220 vs 171–179.** On a crack, |p−1/2| (or |p−1/(n+1)|) is tested at 1e-6, so a slightly wrong p cannot hide. On a *linear notch*, |p−p_Williams| is 1e-6, so Δp=10⁻⁵ FAILs Williams while angular eq is still 10⁻⁵ PASS. On a **power-law notch** there is no p-target at all: only equilibrium/compatibility. That is mechanically correct (p is an eigenvalue) but a 10⁻⁴ error in p on a power-law notch can pass angular eq (tol 10⁻³) if Cartesian also stays under 10⁻³. Measured on the linear 120° notch: Δp=10⁻⁴, Cartesian residual 8.1e-4 PASS if the Williams check were absent.

**Low — no opening / contact.** Negative-KI Williams (closing crack) PASSes every check. Valid for the linear model; wrong for a frictionless contact crack.

**Low — T-stress mixed into Williams at p=1/2.** Mixing ε·T with ε≤10⁻⁴ PASSes (angular residual ~ε). ε=10⁻³ FAILs both equilibria. The leading-term checks do not reject a small regular pollution; they are not meant to.

**Low — F=cos³(θ/2) is not a mutant.** It is proportional to the Williams mode-I Airy function (ratio ≈ 0.75 vs cos(θ/2)+(1/3)cos(3θ/2)). It PASSes because it is correct, not because the checks are blind.

**Handled correctly (not FN/FP):** fake Airy at p=1/2 or HRR p=0.25 n=3; Williams as `law=power n=1` (incompressible linear, J=0.75 PE / 1.0 PS, same as ν=0.5); reported J off by 1.1% FAIL / 0.9% PASS (tol 1%); 150° mode-II linear notch reconstructed from the odd eigenfunction PASSes (Williams mode-I check is n/a).

## 3. Solver constraints / misleading content

**Medium — TOPIC.md:102–113 / crack_tip_checks.py:230–235.** The topic list says checks include J and displacement compatibility; those are implemented for **mode I only**. A mode-II or mixed solver is told the checks will test J; they will not. `reported_J` on a notch is stored but `applies=False`, so a wildly wrong I (1e9) never FAILs.

**Medium — check 7 forces the *leading* J exponent.** A correct higher-order or T-stress field cannot “pass the topic checks” as a package. Fine for Hutchinson’s dominant term; misleading if the next problem is Larsson–Carlsson T-stress or a two-term expansion (TOPIC §5 already names that paper).

**Low — test_checks.py:17–21.** Tests import frozen v1 `hrr` via `PILOTS/corpus/knowledge_v1/...` and a hard-coded `<probe-runs>/P015-build-codex/T2/reference.py`. The checks module itself does not import a solver (good). An agent that treats `test_checks.py` as a template is pulled onto those paths.

**Low — hrr.py** is crack-only, mode I, pure power law, α=π (ENTRY §7 says so). Internal `_` names are marked. Amplitude `alpha` in `amplitude_K` is Ramberg–Osgood α, not the wedge angle. Easy name clash with TOPIC’s α.

**Low — TOPIC.md:36–40.** Power-law strain is written with α, ε₀, σ₀; the checker sets all constants to 1. Stated in the module docstring. A solver that passes dimensional Σ and a dimensional reported J will fail the 1% J match.

**Not a recipe.** TOPIC tells the agent to use any method; ENTRY does not walk through shooting. Normalization “state the convention” is the right constraint.

## 4. Missing for a related problem

For a **linear or power-law notch, mode I**, the topic + checks are almost sufficient (eigenvalue p, faces, eq, compatibility). Missing: mode-II characteristic equation; handling of α=π/2 (no singularity); power-law notch has no independent p-oracle (unavoidable).

For **mode II / mixed mode cracks:** symmetry check exists; J, Iₙ, mixity M^p, and displacement integration do not. Hutchinson 1968 is tensile; Rice–Rosengren / Shih mixity are out of this entry, but TOPIC §5 already lists mode II as in-scope for the equations.

For **T-stress / higher-order terms (P013):** relative Cartesian/Saint-Venant residuals make the valid p=0 field look like equilibrium failure. Need an absolute residual floor (e.g. fail only if |div| > tol × (1 + |grad σ|)).

For **interface cracks (P052/P053):** oscillatory kernels, complex K, energy release rate vs K — not in the checks; scope table only points at the papers.

For **flow theory / growing cracks / finite strain:** scope table is adequate as a warning; no check can replace that.

For **90° corners, contact, or friction:** crash or silent pass (negative KI).

## Verdict

**Usable after the listed fixes.**

The shared mechanics and the Hutchinson entry are sound and usable for a new mode-I crack/notch solver. Do not ship the checks as-is for mode II, T-stress, or a 90° notch: relative residuals and the Williams root-finder will reject or crash correct fields. Minimum fixes: (1) absolute floor on Cartesian and Saint-Venant residuals; (2) return p=0 when the Williams function has no root in (0,1); (3) document that J/displacement checks and “crack exponent fixed by J” apply only to the leading mode-I crack term.
