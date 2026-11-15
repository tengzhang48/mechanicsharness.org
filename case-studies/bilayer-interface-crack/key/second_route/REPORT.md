# Second, independent route for the phase angle of the bilayer interface crack (13 Nov 2026)

**Since the review of draft 0.2 (13 Nov 2026) the key's phase angles are this route's**, set S2 below; the first route
(`../fe_bilayer.py`) is the check. Where this report says "the key", it means the first route, which was the key
when the report was written.

Computed by a separate agent (Claude Opus, a fresh session) that was given the problem and the definitions only. It was
told not to open the key, the scorer or the source paper, and to use another discretisation and another extraction than
the key's. This file is its final report, saved by the organizer (the agent could not write the file itself). The code
in this folder is the agent's. The organizer reran one case (`bilayer.py S2 3 4`, output in
`logs/organizer_rerun_S2_m3.txt`): the interaction integral gives 52.2491° to 52.2492° and 4.2786° to 4.2787° on its
four rings, as below. The plain traction fits on that single mesh scatter more, 52.236° to 52.306° and 4.305° to
4.313°; the "fit" column below is the fit extrapolated over two mesh densities.

## Method

dolfinx 0.10, quadratic (6-node) triangles, plane strain: a polar mesh around the tip (32 m sectors, geometric rings
down to 1e-6 h) blended into a graded triangle grid. The beam stresses are applied as tractions at x1 = −L, the end
x1 = +L is free, three point pins remove rigid motion. K comes from

- (a) a domain interaction integral whose auxiliary field is derived numerically (a complex-exponent eigenproblem for
  the two-material crack tip), averaged over four integration rings; and
- (b) a least-squares fit of the finite-element interface tractions ahead of the tip to
  (σ22 + i σ12) √(2πr) r^(−iε) = K + A1 r + A2 r², extrapolated over two mesh densities.

Neither uses the crack-opening displacement, which is what the key's route uses.

## Results

Units h = 1, E2 = 1, P = 1 (case A), M = 1 (case B), so ψ = arg K. Interaction integral on the finest mesh (m = 4,
truncation L = 4 (h + H)); "fit" is extraction (b). |K| from G uses the beam-theory energy release rate.

| Set | Parameters | ε | ψ_A (°) | ψ_A, fit | ψ_B (°) | ψ_B, fit | \|K_A\| extracted / from G | \|K_B\| extracted / from G |
|---|---|---|---|---|---|---|---|---|
| S1 | E1/E2 = 5, ν1 = 0.22, ν2 = 0.35, h/H = 0.3 | −0.0402107 | 53.954 ± 0.005 | 53.9536 | −19.341 ± 0.005 | −19.3412 | 1.304430 / 1.304430 | 1.563209 / 1.563210 |
| S2 | E1/E2 = 9, ν1 = ν2 = 1/3, h/H = 0.5 | −0.0645318 | 52.249 ± 0.005 | 52.2494 | 4.279 ± 0.005 | 4.2769 | 1.933126 / 1.933127 | 1.629800 / 1.629801 |
| S3 | one material (ν = 0.3), h/H = 0.5 | 0 | 50.944 ± 0.005 | 50.9426 | −22.954 ± 0.005 | −22.9533 | 1.561249 / 1.561249 | 2.598075 / 2.598076 |
| S4 | one material (ν = 0.3), h/H = 0.1 | 0 | 52.043 ± 0.005 | 52.0422 | −37.055 ± 0.005 | −37.0534 | 0.855278 / 0.855278 | 2.450695 / 2.450714 |

Convergence (`logs/summary.md`): mesh densities m = 1 to 4 (3e4 to 7e5 unknowns) and truncations L = 2, 4, 8 times
(h + H) for every set. On the fine meshes each of these changes ψ from the interaction integral by at most 0.0001°:
ring to ring, m = 3 to 4, L = 4 to 8, and another mesh layout. The two extractions differ by at most 0.0018°. ±0.005°
is quoted as a conservative bound. One weak spot: S4, case B, where |K| against G and J against G are off by 1e-5 to
5e-5, believed to be round-off (the thin arm deflects by about 1e4 h); ψ_B is not affected.

## Validation against exact results

1. **One material, h = H.** Case A: ψ_A = 49.10661° against the exact arccos √(3/7) = 49.106605°. Case B: ψ_B =
   0.00000° (pure mode I by symmetry). K_A = 1.732050 + 2.000000 i, which is √3 + 2i.
2. **Finite interface crack between two half-planes** (the S2 materials, ε = −0.0645; a plate with the remote
   tractions on its boundary; K = (σ + iτ)(1 + 2iε) √(πa) (2a)^(−iε) at the right tip). Phase error of the interaction
   integral at most 0.001° for every plate size and m ≥ 2. Error of |K|: 1.9e-3, 4.7e-4, 1.2e-4 for plate half-widths
   25, 50, 100 crack half-lengths, i.e. proportional to (a/W)², the finite-plate effect. With the S1 materials and with
   them exchanged (ε = +0.040) the phase error is at most 0.0001° (`logs/griffith_table.md`).
   The convention was checked two ways: the expansion of the closed-form interface traction has the r^(+iε) form at
   the right tip; and numerically the opposite sign of ε would put the phase at +4.79° instead of −4.79°.
3. **Relative phase from energy alone.** K is linear in the loads and G = a P² + 2 b P M + c M², so
   cos(ψ_A − ψ_B) = b/√(ac) exactly. The finite-element phases satisfy this to 0.0000° for every set: ψ_A − ψ_B =
   73.2947° (S1), 47.9705° (S2), 73.8979° (S3), 89.0979° (S4) (`logs/energy_phase_check.md`). This fixes the relative
   sign of K2 between the two load cases, not the absolute phase.

## Conventions

x1 ahead of the tip, x2 into layer 1, material 1 above. K is defined only by σ22 + i σ12 = K r^(iε)/√(2πr) on the
interface ahead of the tip; r^(iε) = exp(iε ln r) with ln r real; ψ = arg(K h^(iε)) with h = 1, in (−180°, 180°].
For the eigenvalue 1/2 + iε, with ε from the task's formula, the eigenfield has σ_rθ/σ_θθ = −i on the interface for
both signs of ε, so the task's formula for ε and "material 1 above" are consistent. Physical sign for one material
and h = H: behind the tip the upper face slides towards +x1 relative to the lower face, so K2 > 0.

## What this route found about the task

With the S2 parameters both ψ_A = 52.25° and ψ_B = +4.28° are positive, so under positive P and M the phase angle
stays between them and never reaches zero: the question "the ratio M/(P h) at which ψ = 0" had no answer. (Im K = 0
only at M/(P h) = −12.57, where Re K < 0.) The question was replaced by the phase angle under both loads with
M = P h; from this route's two K it is 30.433°.

## Files

`asym.py` (the crack-tip eigenfield and auxiliary fields, with self-tests), `check_aux_J.py`, `check_eigs.py`,
`meshgen.py`, `fe.py` (the solve and both extractions), `bilayer.py` (one run: `python bilayer.py S2 3 4`),
`run_all.sh`, `summarize.py` (builds `logs/summary.md` and `results/final.json`), `griffith.py` (validation 2),
`energy_phase_check.py`, `robust_R0.py`. Run with the Python of the dolfinx environment (`fenicsx-env`); the largest
run (1.1 million unknowns) takes about a minute.
