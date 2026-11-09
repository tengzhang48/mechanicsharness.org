# Task: build three short probe problems from a mechanics paper

Paper: Hutchinson J. (1968). Singular behaviour at the end of a tensile crack in a hardening material. Journal of the Mechanics and Physics of Solids 16, 13-31 (corpus id P015). Page images `pages/page-NN.png` (300 dpi) and the text
layer `paper.txt` are in this folder.

Another AI agent (from a different model family) will later try to solve your problems, once with a
knowledge unit about this paper and once with the problem statement only. Your problems test
whether knowledge of this paper **transfers to related work**, not whether the agent can recite it.
You will not see the knowledge unit, and the solver will not see the paper.

Write three probes, each in its own folder `T1/`, `T2/`, `T3/`, each with `TASK_PROBLEM.md`,
`REFERENCE.md` (with your reference code in the same folder) and `RUBRIC.md` (100 points):

- **T1, parameter transfer.** Ask for a quantity of the paper's model at inputs the paper does not
  tabulate or state (for example another parameter value or another regime of the same model).
  The answer must be a number (or a few numbers) with a stated tolerance.
- **T2, variant problem.** A neighbouring problem that needs the paper's mechanics *adapted*: a
  different geometry, loading, boundary condition or material arrangement, in the spirit of "a
  related new problem". It must be solvable by a capable agent in **45 minutes** with Python
  (numpy, scipy). Ask for numbers plus a short justification of the model chosen.
- **T3, applicability.** Three short cases: one where the paper's model applies, one where it does
  not, and one borderline. For each, the solver must decide and give the physical reason. A good
  answer may adapt the model, name missing inputs, or propose additional physics. **Penalize
  unnecessary refusal as much as wrongly applying the model.**

Requirements:
- **Self-contained statements:** give every input, convention and unit the solver needs. Do not
  mention the paper's title, authors or equation numbers, and do not cue the answer or the model
  choice (for T2 and T3 the choice is part of what is tested).
- **Independently checked references:** compute each reference with your own code and check it by
  a second route (a limit case, another method, or a known special solution). Report both in
  `REFERENCE.md` with the agreement achieved. A reference you could not check must say so.
- **Rubrics** reward correct model choice, correct adaptation, verification and scope recognition;
  code reuse is not required. State partial credit explicitly.
- Keep each problem short enough for a 45-minute budget.

Also write `BUILDER_NOTES.md` (what each probe tests, how it relates to the paper, anything a
grader should know) and keep `PROCESS.md` (one line per step, each starting with a time pasted
from `date -u +%H:%M:%S`; never type a time). Rules: work only inside this folder (use
`./scratch/`); no web search or fetch. Time limit: 90 minutes.
