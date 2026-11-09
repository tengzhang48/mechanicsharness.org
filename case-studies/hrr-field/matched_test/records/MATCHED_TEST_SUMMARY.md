# Matched test of the entry format: results (GLM and DeepSeek)

2026-10-07, logs P251–P260. Reproduce the tables with `python3 analyze_matched.py`
(output saved in `matched_analysis_output.txt`; per-run data in `matched_results.json`).
**Two runs per cell; descriptive only.**

## Design

- Knowledge: the frozen crack-tip asymptotics topic v0.2 with the Hutchinson 1968 entry
  (`library/topics/crack_tip_asymptotics`), frozen before the probes were read (P249–P250).
- Conditions: **B0** topic + checks + entry, no solver code; **B1** the same plus the entry's
  reference code (only solver access differs: 10 lines of the entry and one file); **C** nothing.
- Probes, built by Codex from the paper alone, references checked by two routes (V1 also by Opus,
  independently): **V1** power-law antiplane (mode III) crack, outside the topic's in-plane scope;
  **V2** inverse calibration of a bilinear material in plane stress; **A1** three applicability cases.
- Solvers: GLM-5.3 and DeepSeek-V4.1-Flash (opencode), 2 repeats: 36 runs, all ended normally with
  an answer, all audited (no web use, no reference files; FAILURES 58: one run wrote into another
  run's folder and tried to remove it, harmlessly). Graded blind by Codex; boundary checks coded
  separately and blind to the grade (`BOUNDARY_CODING.md`, pre-registered).
- Kimi K3 (the stronger solver) could not run in parallel (Azure token-rate limit, P253–P255); after
  the ceiling was seen it was stopped (P258). Its one completed cell is not analysed.

## Results

| Probe | Solver | B0 | B1 | C |
|---|---|---|---|---|
| V1 | GLM | 100, 100 | 99, 100 | 100, 100 |
| V1 | DeepSeek | 100, 100 | 100, 100 | 100, **76** |
| V2 | GLM | 98, 100 | 100, 100 | 100, 100 |
| V2 | DeepSeek | 100, 100 | 100, 100 | 100, 100 |
| A1 | GLM | 93, 99 | 93, 96 | 91, 91 |
| A1 | DeepSeek | 93, 93 | 93, 93 | 94, 94 |

Means: B0 98.0, B1 97.8, C 95.5. Run time: GLM about 14–15 min in every condition; DeepSeek
3–5 min (knowledge adds about 1.5 min of reading).

1. **The probes sit at the ceiling** for both solvers in every condition. Fully specified
   variants of a classical crack-tip problem are within the reach of these agents without any
   knowledge. This is a result in itself, and the reason the next round moves to research-level
   problems from development papers (P258).
2. **The one real failure was an execution error in J**, not a missing boundary check: DeepSeek
   without knowledge computed the J integral with a derivative-parity error, so the amplitude and
   every downstream number were off by a consistent factor (grader M20). With the topic or entry,
   DeepSeek scored 100 in all four V1 runs; the topic warns about J computation and its checks
   compute J independently. The same error class appeared in the probe pilot (HRR T1) and in the
   HRR unit itself.
3. **Applicability (A1)**: a small gradient for GLM (B0 96, B1 94.5, C 91), none for DeepSeek.
4. **Use of the knowledge**: the checks were imported by the solver's own code in 9 of 24
   knowledge runs (B0 5, B1 4); the reference code in 3 of 12 B1 runs. On V1 (mode III, outside
   the topic's in-plane scope) the solvers derived their own equations.

## Is the boundary check the key? (pre-registered question)

| Condition | Existence (B1) | Limits (B2) | Applicability (B3) | Degeneracy (B4) | ... of which done **before** the numerics (B1/B2/B3/B4) |
|---|---|---|---|---|---|
| B0 | 11/12 | 8/12 | 8/12 | 1/12 | 2 / 2 / 0 / 0 |
| B1 | 12/12 | 9/12 | 9/12 | 2/12 | 4 / 2 / 0 / 1 |
| C | 9/12 | 6/12 | 10/12 | 0/12 | 3 / 0 / 0 / 0 |

- **Against the pre-registered "support" criterion, this round gives no support:** the only
  failure had correct existence, limit and applicability checks; its error was in execution.
  With scores at the ceiling the test is weak either way.
- **What the coding does show**: with knowledge, existence and limit checks are somewhat more
  frequent (B1 condition: 12/12 and 9/12 vs 9/12 and 6/12 without), and degeneracy is considered
  only with knowledge (rarely). **Most boundary checks are done after the numerics, as final
  verification, not before**, in every condition; the knowledge shifts this only slightly. Agents
  do not, on their own, follow the "boundary first" practice (P245).
- Caveat: B1/B2 coding is partly not applicable to A1 (an applicability task), so per-probe
  differences confound the score-vs-checks comparison.

## Limits

Two runs per cell; ceiling scores; partial blinding (some answers mention the knowledge);
Codex both built the probes and graded; one topic, one paper; the K3 arm missing.

## Next

Research-level tasks from development papers in the nonlinear crack-tip area (P088 first),
calibrated in condition C so that they are not at the ceiling (P258). Whether a "boundary
first" instruction changes behaviour is a separate, testable intervention.
