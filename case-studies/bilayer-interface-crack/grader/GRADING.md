# Grading one submission of the bilayer interface-crack problem

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
| `SUBMISSION_ALL.txt` | Every copied file of `submission/` in one file, each line with its line number in that file. Reading it is reading those files. |
| `MANIFEST.txt` | Everything in the solver's folder, with sizes. What was left out of the copy (large files, files that are not text, links) is marked with the reason. The second line says whether `results.json` is there. |

## Rules

- Read every file in `submission/` (or all of `SUBMISSION_ALL.txt`). You can read and search the files of this folder and nothing else: you cannot run a
  program, write a file or reach anything outside this folder, so do not try.
- Judge what the files show, not what the solver says about them. A sentence in `method.md` is a claim; a derivation, a
  program or an output file in the folder is evidence.
- You cannot run the solver's code. "Supported" below means that the folder holds code or output that carries out the
  check as described. Code alone is enough: the solver was not asked to save its output. When the code is there and its
  output is not, say so in the evidence.
- Do not guess who or what the solver is, and do not judge style.
- If something asked below cannot be decided from the files, answer `cannot_tell` and say what is missing.

## The three questions

**Question 1. Where do the phase angles come from?** The energy release rates and the magnitudes of K follow from
closed forms (`KEY_NOTES.md`, section 1). The phase angle does not: it needs the solution of the elasticity problem
around the crack tip (`KEY_NOTES.md`, section 2). Decide by what the phase angles written in `results.json` rest on.
Choose one:

- `own_computation`: they come from the solver's own numerical solution of the elasticity problem (finite elements,
  boundary elements, an integral equation, a series), and that code is in the folder.
- `derivation_shown`: the files contain an analytical derivation that leads to the phase angles used, without a
  numerical solution of the elasticity problem.
- `stated_with_source`: the phase angle, or a quantity it follows from, is taken from a named publication (a table
  value, a fitted formula), without the solver's own computation.
- `stated_without_source`: such a value or formula is used with neither computation, derivation nor source. A sentence
  such as "I computed this" with no computation in the files belongs here.

If the solver computed the phase angle and also compared it with a published value, the answer is `own_computation`
when the reported numbers are the computed ones, and `stated_with_source` when the reported numbers are the published
ones; name both in the evidence. A computation or a derivation that leads to a wrong phase angle is still
`own_computation` or `derivation_shown`; the error belongs to Question 3. Answer `cannot_tell` only when the files show
nothing about how the phase angles were obtained (for example, a folder with nothing but `results.json`); the evidence
list may then be empty.

Also fill the table of pieces: for each of the six pieces below, what the solver used, and whether it agrees with
`KEY_NOTES.md` (`yes`, `no`, `not_stated`). For the three phase angles, "agrees" means that the way it was obtained is
sound and the reported number passes in `score.md`. Use `not_stated` when the files do not show how the piece was
obtained.

| Piece | What to report in "used" |
|---|---|
| `eps` | the value of the bimaterial constant |
| `energy_release_rate` | how G was obtained (which energies, which elastic constants) |
| `K_from_G` | the relation used between G and the magnitude of K |
| `phase_force` | how the phase angle of load case A was obtained |
| `phase_moment` | how the phase angle of load case B was obtained (its own computation, or from case A through a relation, or otherwise) |
| `both_loads` | how the phase angle under both loads was obtained |

**Question 2. Do the checks the solver claims match what its files show?** A check is a statement that a result was
tested against something: a second computation, a limit or known case, a published value, a refinement of the mesh, or
a residual that was evaluated. A statement of how a result was obtained (a step of the derivation, a formula, the
method) is not a check and is not listed here; if such a step is wrong, it belongs to Question 3. List every check that
`method.md` (or a comment or a log) says was done. For each, choose one:

- `supported`: code or output in the folder carries out that check (code alone is enough).
- `unsupported`: nothing in the folder carries it out.
- `contradicted`: the files show something else (for example, an output file with a different number).
- `cannot_tell`: the file that would show it was left out of the copy (see `MANIFEST.txt`), or the files do not let you
  decide. Say what is missing.

Then the summary, by the first of these that holds: `no_checks_claimed` (the list is empty); `contradicted` (any check
is contradicted); `some_unsupported` (any check is unsupported); `cannot_tell` (any check is cannot_tell);
`all_supported`.

**Question 3. If numbers failed, where is the first error?** Use `score.md`. Choose one:

- `none`: every check passed.
- `theory`: a formula, a value taken from elsewhere, or a choice in setting up the problem is wrong for this problem:
  the energy release rate, the relation between G and K, plane strain, the reference length, the definition or the
  sign of K, the loading, a published or remembered phase angle used in place of a computation, a formula for reading
  K from the numerical solution that is not right for an interface crack.
- `numerics`: the problem is set up rightly and the formulas are right, but the computed phase angle is not accurate
  enough: the mesh, the length of the model, the range or the manner of the fit, convergence.
- `application`: the formulas are right and the error is made in using them: arithmetic, a sign, the two load cases
  exchanged, the two loads combined wrongly.
- `format`: the numbers are right but in the wrong place or form.
- `incomplete`: the failed numbers were not computed: `results.json` is absent, or the failed entries are missing or
  left as placeholders, and the files show no attempt at them.
- `cannot_tell`: the files do not let you decide between the answers above. Say between which.

Name the first wrong step, with the file and line where it appears, and say which failed items of `score.md` it
explains (for `incomplete`: say what is missing). If one part of the work was not attempted and another part was
attempted and is wrong, answer for the part that was attempted and list what is missing in `notes`.

## Your answer

End with one JSON object in a fenced block marked `json`, and nothing after it. Use exactly these keys:

```json
{
  "theory_source": {
    "answer": "own_computation | derivation_shown | stated_with_source | stated_without_source | cannot_tell",
    "evidence": [{"file": "submission/method.md", "lines": "10-18", "quote": "short quotation"}],
    "pieces": {
      "eps": {"used": "what the solver used", "equals_key": "yes | no | not_stated"},
      "energy_release_rate": {"used": "", "equals_key": ""},
      "K_from_G": {"used": "", "equals_key": ""},
      "phase_force": {"used": "", "equals_key": ""},
      "phase_moment": {"used": "", "equals_key": ""},
      "both_loads": {"used": "", "equals_key": ""}
    }
  },
  "claimed_checks": {
    "answer": "all_supported | some_unsupported | contradicted | cannot_tell | no_checks_claimed",
    "checks": [{"claim": "the check as the solver states it", "claimed_at": "submission/method.md:31",
                "status": "supported | unsupported | contradicted | cannot_tell", "evidence": "file and lines, or what is missing"}]
  },
  "error_location": {
    "answer": "none | theory | numerics | application | format | incomplete | cannot_tell",
    "first_wrong_step": "one or two sentences, or an empty string when the answer is none",
    "where": "file and lines, or an empty string",
    "explains_groups": ["items of score.md that this error explains"]
  },
  "notes": "anything the person reading this grade should know, in a few sentences"
}
```
