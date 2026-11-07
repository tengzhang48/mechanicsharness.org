# Decision: disclose coordinator text without changing the statement

**Recommendation / inference:** approve a dated disclosure in FIXES.md and on the private case page. Preserve the exact
statement as run. A clean solver statement would be a separately versioned task with new runs. This is a metadata
decision, not a key correction, and none of those changes has been applied in this step.

**Read:** `task/PROBLEM.md:3–4` says:

> Draft for Teng's review. This is the proposed source-faithful reproduction
> task; source solutions and the proposed geometry extension are excluded.

**Read:** `FIXES.md:6–13,35–37` already names the defect and preserves it. The statement also leaves numerical accuracy
unspecified; `FIXES.md:21–32` records the later 0.5% checker alongside the historical 10⁻⁶ count.
`audit_2026-10-25/AUDIT.md`, §1, treats the coordinator wording as part of the input, not a post-run addition.

**Inference:** the cue tells solvers the task reproduces a source, and its “proposed” wording makes the input look less
final. It does not by itself prove that a solver recalled the paper or establish a numerical score adjustment. Removing
it now would misrepresent what was tested. Keep this disclosure adjacent to run results, where readers can see it before
interpreting the results.

**Proposed FIXES.md addition / inference** (append under “Statement: kept, flaws named”):

> **Release annotation, 29 October 2026.** The solver input included the coordinator's draft/review wording at
> task/PROBLEM.md lines 3–4, quoted above. Those lines were part of the statement as run and remain unchanged. Treat these
> runs as source-reproduction trials with that cue; they are not evidence from a clean, uncued statement. This annotation
> changes neither the supplied problem nor the historical grades.

**Proposed case-page note / inference** (next to the task link and before results):

> **Statement as run.** Agents received two coordinator lines identifying the task as a draft source-faithful
> reproduction. The archived task preserves those lines. A cleaned statement would require a new version and new runs;
> these results refer to the original input. Numerical agreement under the original 10⁻⁶ comparison and the later 0.5%
> checker is reported separately.

**Read / source location update:** the corpus filename search found no Xia–Hutchinson paper. Following the local path in
`audit_2026-10-25/AUDIT.md:5–6,8` did locate
`<source paper PDF, not included>`.
Its first rendered page verifies the title/authors for the outreach ledger. This packet needs no source equation to
decide the input-preservation issue and does not claim a fresh source verification of the film solution.

**Exactly if approved:** append the above annotation to `FIXES.md`; add the case note to the embedded Xia–Hutchinson
case page in `case_library/AI_Agents_in_Fracture_Mechanics_2026-10-26.html` (and to its case source if/when that source
is restored). Link this decision from the case note. Do not alter `task/PROBLEM.md`, the copies in run workspaces, any
key/rubric/checker or the old scores. If an uncued task is wanted, create `xia_hutchinson_film_v2/` with a cleaned task
and explicit accuracy requirement, freeze it and run it anew.
