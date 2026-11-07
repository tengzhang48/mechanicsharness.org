# Grading process

19:01:23 UTC — Read TASK.md and inventoried the submission. Read the problem, corrected key, automatic checks, and initial submission files. Work stays in this folder; no web access. Will inspect all claims and derivations, reproduce code on a scratch copy, and write the required grading profile.

19:01:51 UTC — Inspected the full answer, analytic script, verification scripts, results, and saved transient output. Identified seven matching values, the key's explicit partial-credit rule for K1_0 = 1, an incorrect short-time coefficient, and a written-versus-saved transient-K mismatch. Prepared a scratch copy for reproduction. The analytic and static-verification reruns succeeded; transient rerun started (tool session 65040).

19:02:22 UTC — Confirmed solve.py rewrites results.json identically to the submission and reproduces all nine printed values. Inspected the original NPZ archive: K values are negative correction-field values, inconsistent with the answer's quoted positive total values. Computed deviations from the corrected key: K1_0 +23.82%, short-time coefficient -53.76%.

19:04:18 UTC — Completed independent diagnostics in scratch/audit.py. Demonstrated an incorrect volume integral for nonuniform loading (2.932153 coded versus 2.513274 integrated), missing L in the stress-intensity quadrature, malformed implicit Neumann rows, a one-cell reflection offset, and incorrect off-crack boundary data in the steady strip pressure. Verified the conditional diffusivity algebra. Transient rerun finished with the same V, K, and A_fit as the saved output.

19:06:42 UTC — Added a grid/box diagnostic: the steady FFT trace error persists as resolution increases and falls as the box expands, contradicting attribution to tip aliasing. Classified submitted source attributions, six questions/assumptions, seven reported checks, and the unsupported 3D-trace interpretation challenge.

19:08:56 UTC — Wrote GRADE.md with all nine requested sections, corrected-key answer judgments (a partly, b right, c partly), requested supporting-state credit, three decisive errors, ten side errors, two presentation errors, reproduction findings, and the required final JSON block. Included the key's initial-K order-of-limits distinction and a dimensionally incorrect prefactor in its historical opening–K transcription; neither changes the final question-level profile.

19:09:54 UTC — Final validation passed: nine report sections, parseable final JSON with nothing following it, matching numbered error counts, exact results.json reproduction, unchanged original submission files, and reproduced transient arrays within 1e-12. GRADE.md and the process log are complete; all scratch evidence remains inside this folder. No external sources or processes were accessed or stopped.
