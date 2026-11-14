# Grading one submission of the trapped-front problem

You are grading one solver's work on a fracture mechanics problem. The numbers have already been checked by a program;
its report is `score.md`. Your job is to read the solver's files and answer three questions about **how** the numbers
were obtained. A person will read your grade and make the final judgment, so every answer must point to the lines it
rests on.

## What is in this folder

| Path | What it is |
|---|---|
| `TASK.md` | The problem exactly as the solver received it. The solver received nothing else. |
| `KEY_NOTES.md` | The right theory and the usual mistakes. The solver never saw this. |
| `score.md` | The program's check of the solver's numbers against the reference values. |
| `submission/` | A copy of the solver's folder: `results.json`, `method.md`, and any code and output files that are text. |
| `MANIFEST.txt` | Everything in the solver's folder, with sizes. What was left out of the copy (large files, files that are not text, links) is marked with the reason. The second line says whether `results.json` is there. |

## Rules

- Read every file in `submission/`. You can read and search the files of this folder and nothing else: you cannot run a
  program, write a file or reach anything outside this folder, so do not try.
- Judge what the files show, not what the solver says about them. A sentence in `method.md` is a claim; a derivation, a
  program or an output file in the folder is evidence.
- You cannot run the solver's code. "Supported" below means that the folder holds code or output that carries out the
  check as described. Code alone is enough: the solver was not asked to save its output. When the code is there and its
  output is not, say so in the evidence.
- Do not guess who or what the solver is, and do not judge style.
- If something asked below cannot be decided from the files, answer `cannot_tell` and say what is missing.

## The three questions

**Question 1. Where does the first-order theory come from?** The problem needs the first-order change of the three
stress intensity factors along a slightly curved front: a nonlocal term for each mode and a coupling of modes 2 and 3
through the slope of the front (`KEY_NOTES.md`, section 1). Choose one:

- `derivation_shown`: the files contain a derivation from the equations of elasticity that leads to the coefficients used.
- `own_computation`: the first-order response (the coefficients, or the stress intensity factors along the front
  directly) comes from the solver's own numerical solution of the three-dimensional crack problem, and that code is in
  the folder.
- `stated_with_source`: the formulas are stated and attributed to a named publication, without derivation.
- `stated_without_source`: the formulas are stated with neither derivation nor source. A sentence such as "I derived
  this" with no derivation in the files belongs here.

If more than one applies, choose the first in the order above: a derivation that is shown counts even when a source is
also cited (name the source in the evidence). A derivation or a computation that leads to wrong coefficients is still
`derivation_shown` or `own_computation`; the error belongs to Question 3. Answer `cannot_tell` only when the files show
no formulas at all (for example, a folder with nothing but `results.json`); the evidence list may then be empty.

Also fill the coefficient table: for each of the six coefficients of `KEY_NOTES.md`, section 1, the value or expression
the solver used and whether it equals the key's (`yes`, `no`, `not_stated`). Use `not_stated` when the solver never
forms that coefficient (for example, in a direct numerical solution).

**Question 2. Do the checks the solver claims match what its files show?** A check is a statement that a result was
tested against something: a second computation, a limit or known case, a published value, or a residual that was
evaluated. A statement of how a result was obtained (a step of the derivation, a formula, the method of a search) is
not a check and is not listed here; if such a step is wrong, it belongs to Question 3. List every check that
`method.md` (or a comment or a log) says was done. For each, choose one:

- `supported`: code or output in the folder carries out that check (code alone is enough).
- `unsupported`: nothing in the folder carries it out.
- `contradicted`: the files show something else (for example, an output file with a different number).
- `cannot_tell`: the file that would show it was left out of the copy (see `MANIFEST.txt`), or the files do not let you
  decide. Say what is missing.

Then the summary, by the first of these that holds: `no_checks_claimed` (the list is empty); `contradicted` (any check
is contradicted); `some_unsupported` (any check is unsupported); `cannot_tell` (any check is cannot_tell);
`all_supported`. The summary is one of these five words; it is not a status word.

**Question 3. If numbers failed, where is the first error?** Use `score.md`. Choose one:

- `none`: every check passed.
- `theory`: the first-order theory itself is wrong (a coefficient of `KEY_NOTES.md`, section 1, or the energy release rate).
- `application`: the theory is right, and the error is made in using it (the operator on a cosine, the Fourier
  inversion, a sign, a normalization, the sampling points, the search for the maximum).
- `format`: the numbers are right but in the wrong place or form.
- `incomplete`: the failed numbers were not computed: `results.json` is absent, or the failed entries are missing or
  left as placeholders, and the files show no attempt at them.
- `cannot_tell`.

Name the first wrong step, with the file and line where it appears, and say which failed groups of `score.md` it
explains (for `incomplete`: say what is missing). If one part of the work was not attempted and another part was
attempted and is wrong, answer for the part that was attempted and list what is missing in `notes`.

## Your answer

End with one JSON object in a fenced block marked `json`, and nothing after it. Use exactly these keys:

```json
{
  "theory_source": {
    "answer": "derivation_shown | own_computation | stated_with_source | stated_without_source | cannot_tell",
    "evidence": [{"file": "submission/method.md", "lines": "10-18", "quote": "short quotation"}],
    "coefficients": {
      "p1": {"used": "what the solver used", "equals_key": "yes | no | not_stated"},
      "p2": {"used": "", "equals_key": ""},
      "p3": {"used": "", "equals_key": ""},
      "b2": {"used": "", "equals_key": ""},
      "b3": {"used": "", "equals_key": ""},
      "m3": {"used": "", "equals_key": ""}
    }
  },
  "claimed_checks": {
    "answer": "all_supported | some_unsupported | contradicted | cannot_tell | no_checks_claimed",
    "checks": [{"claim": "the check as the solver states it", "claimed_at": "submission/method.md:31",
                "status": "supported | unsupported | contradicted | cannot_tell", "evidence": "file and lines, or what is missing"}]
  },
  "error_location": {
    "answer": "none | theory | application | format | incomplete | cannot_tell",
    "first_wrong_step": "one or two sentences, or an empty string when the answer is none",
    "where": "file and lines, or an empty string",
    "explains_groups": ["group titles of score.md that this error explains"]
  },
  "notes": "anything the person reading this grade should know, in a few sentences"
}
```
