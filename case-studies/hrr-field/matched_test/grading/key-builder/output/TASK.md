# Task: build three fresh probe problems from a mechanics paper

Paper: Hutchinson J. (1968). Singular behaviour at the end of a tensile crack in a hardening material. Journal of the Mechanics and Physics of Solids 16, 13-31 (corpus id P015). Page images `pages/page-NN.png` (300 dpi) and the text
layer `paper.txt` are in this folder.

Other AI agents (from other model families) will later try to solve your problems under
different conditions: with knowledge material about this paper's subject, or with the problem
statement only. Your problems test whether mechanics knowledge **transfers to related new
problems**: choosing the governing mechanism, formulating and adapting a model, solving and
verifying it, and recognizing when the available physics is inadequate. You will not see the
knowledge material, and the solvers will not see the paper.

**These must be fresh.** Earlier probes from this paper (already used; do not repeat or lightly
vary them): a plane-strain crack-tip field at hardening exponent n = 6.5 (stress exponent, a
near-tip stress amplitude H, and an energy ratio); a 120-degree notch in a plane-strain n = 4
material (stress exponent and a stress ratio at the symmetry line); an applicability set about
cyclic loading, large-scale yielding and J-contour agreement.

Write three probes, each in its own folder `V1/`, `V2/`, `A1/`, each with `TASK_PROBLEM.md`,
`REFERENCE.md` (with your reference code in the same folder) and `RUBRIC.md` (100 points):

- **V1 and V2, variant problems.** Natural neighbouring problems a researcher would pose that
  need the paper's mechanics *adapted* (another geometry, loading mode, constraint, material
  law or boundary condition). Each must be solvable by a capable agent in **45 minutes** with
  Python (numpy, scipy). Ask for numbers with stated tolerances plus a short justification of
  the model chosen. At least one of V1/V2 should be a case where the paper's model is useful
  once adapted.
- **A1, applicability.** Three short cases: one where the model applies (possibly after
  adaptation), one where it does not and additional physics is needed (say which), and one
  borderline. For each the solver decides and gives the physical reason. Penalize unnecessary
  refusal as much as wrongly applying the model.

Requirements:
- **Self-contained statements** with every input, convention and unit. Do not mention the
  paper's title, authors or equation numbers, and do not cue the answer or the model choice.
- **Independently checked references:** compute each reference with your own code and check it
  by a second route (a limit case, another method, a known special solution); report both and
  the agreement in `REFERENCE.md`. A reference you could not check must say so.
- **Rubrics** score model choice and its justification, numerical correctness, completion,
  verification, and recognition of missing physics, with explicit partial credit. Code reuse is
  not required; any correct method earns full credit.

Also write `BUILDER_NOTES.md` (what each probe tests, how it relates to the paper, anything a
grader should know) and keep `PROCESS.md` (one line per step, each starting with a time pasted
from `date -u +%H:%M:%S`; never type a time). Rules: work only inside this folder (use
`./scratch/`); no web search or fetch. Time limit: 90 minutes.

Never stop or kill a process you did not start. Kill only by the PID you recorded when you started it; never use pkill, killall, or pgrep-based loops (this machine is shared with other runs).
