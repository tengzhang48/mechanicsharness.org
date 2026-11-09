# Transfer-probe pilot: does a frozen knowledge unit help an agent on a related problem?

> **Excerpt for the HRR field case, made on 8 November 2026.** Kept, unchanged: the title, the design and the rows of the HRR problems (P015) in the score tables, with their headers. Left out: the rows and findings of the other problems of the pilot (P059, P038) and the later sections.

| Paper | Probe | GLM | DeepSeek-Azure |
|---|---|---|---|
| P015 | T1 parameter | 97 / **0** (timeout, no answer) | 100 / 55 |
| P015 | T2 variant | **0** (timeout, no answer) / 98 | 88 / 76 |
| P015 | T3 applicability | 94 / 88 | 83 / 70 |

| Paper | Mean with | Mean without |
|---|---|---|
| P015 (HRR) | 77.0 | 64.5 |

## Design

- **Units:** three frozen knowledge units (`corpus/knowledge_v1/`): Evans & Hutchinson 1984
  (P059, thin-film buckling delamination), Bao & Suo 1992 (P038, crack bridging),
  Hutchinson 1968 (P015, the HRR field; frozen with listed caveats after four correction rounds).
- **Probes per paper**, built by Codex from the paper alone (it never saw the unit), each with a
  reference, a rubric and an independent second route:
  T1 = the paper's method with new parameters; T2 = a variant the paper does not solve;
  T3 = applicability judgement (one valid, one invalid, one borderline case), where an
  unnecessary refusal is capped like a wrong application.
- **Conditions:** with the unit in `./knowledge/` vs without. **Solvers:** GLM-5.3 (opencode) and
  DeepSeek-V4.1-Flash (Azure, Cline). One run per cell; 50 min limit; detached; every session
  audited (no web use, no run touched the reference files).
- **Grading:** fresh Codex sessions, blind labels (`LABELS_coordinator_only.json`), rubric with
  partial credit; graders rerun the solver's code where numbers depend on it. A run with no
  ANSWER.md at the time limit scores 0 and is not graded.
