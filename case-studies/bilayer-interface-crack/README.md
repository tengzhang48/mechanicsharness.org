# Interface crack between two elastic layers: the mode mix

A MechanicsHarness case. The case page explains the mechanics, the solution and what each agent did:
https://mechanicsharness.org/case-studies/bilayer-interface-crack.html

The problem comes from Z. Suo and J. W. Hutchinson (1990), "Interface crack between two elastic layers", *Int. J. Fract.*
43, 1–18: a crack on the interface of two bonded layers, loaded far behind the tip by a force and a moment in the upper
layer. The parameters (E₁ = 9E₂, ν₁ = ν₂ = 1/3, h/H = 0.5) are our choice. The agents received the problem without the
paper. The paper is not included here.

This folder lets you attempt the task in either version, check a new answer, and check every saved answer and grade
against the key.

## The two versions

| Version | Task | What it gives | Runs |
|---|---|---|---|
| With the definitions | `task/PUBLIC_TASK.md` | the formula for the bimaterial constant ε, the definition of the complex stress intensity factor K with its factor r^(iε), and the phase angle ψ with the reference length h | `runs/*_1` |
| Without them | `task/VARIANT_NO_DEFINITIONS.md` | none of these; it asks the amplitude of the crack-tip traction, its phase at r = h and how that phase changes when r grows tenfold | `runs/*_nodef1` |

Sections 1 and 2 (the layers and the loads) are word for word the same. `python3 make_workspace.py WS no_definitions`
builds a workspace of the second version.

## The answer in one paragraph

G follows from the beam energies far behind the tip: G E₂h/P² = 287/162 under the force, G E₂h³/M² = 34/27 under the
moment. ε = ln(2/3)/2π = −0.0645, so the phase of the crack-tip traction falls by 8.514° for every tenfold increase of the
distance from the tip. |K| follows from G with the two-material relation (cosh²(πε) = 25/24). The phase angle at r = h has
no closed form: from two independent finite-element computations, 52.249° under the force and 4.279° under the moment
(the two agree within 0.02°), and 30.433° under both with M = Ph. The energy alone fixes the size of their
difference, 47.97°. At this point the paper's printed table gives 51.4°; the difference is not explained (the first
computation against the printed tables: `key/validate_fe.log`; for one material, against the authors' later formula:
`key/check_homogeneous.log`).

## What is here

| Path | Content |
|---|---|
| `task/` | The two versions of the task, exactly as given. |
| `key/reference.json`, `key/reference_bilayer.py` | The reference of the first version: closed forms, and the phase angles from the second computation, with the first computation's values and their differences recorded beside them. |
| `key/reference_no_definitions.py`, `.json` | The same reference in the terms of the second version. |
| `key/second_route/` | The second finite-element computation (quadratic triangles, dolfinx), by a separate agent that saw neither the first computation nor the paper: its report, code, results and convergence tables. |
| `key/fe_bilayer.py`, `key/validate_fe.py`, `key/check_*.py` and their logs | The first computation (NumPy, SciPy) and its checks against exact cases and the paper's printed tables. |
| `key/container_timing.py`, `.log` | The first computation run in the agents' container. |
| `scorer/` | The eight checks of each version, the points by piece of physics (approved 14 November 2026), and the controls of each version. |
| `grader/` | The reading of each run's files by Grok 4.7: the instructions of each version, its form, and two control folders. The tools themselves are those of the trapped-front case. |
| `runs/<run>/` | Each run: the facts of the run, `AUDIT.md`, `submission/` (every file it saved up to 2 MB, as delivered), `score/`, `grade/`, and `messages.md` (its visible messages and tool steps). |
| `runs/RESULTS.md` | The results of both versions and what they show. |
| `runs/launch_*.sh`, `runs/process_run.sh` | How the runs were launched and processed. |
| `check_package.py` | Checks this folder (below). |
| `RELEASE_EDITS.md`, `SHA256SUMS`, `LICENSE`, `LICENSE-CONTENT` | Changes made for the release, hashes, licences. |

Not included: the full event streams of the runs and of the grader, the source paper, and the organizer's private notes.

## The runs

| Run | Model | Version | Time | Checks | Points |
|---|---|---|---|---:|---:|
| `codex_1` | GPT-6.1 Sol, Codex CLI, reasoning max | with the definitions | 11.6 min | 8 | 100 |
| `kimi_1` | Kimi K3.1 Pro, opencode, max | with the definitions | 8.5 min | 8 | 100 |
| `deepseek_1` | DeepSeek V4.1 Flash, opencode, max | with the definitions | 41.6 min | 8 | 100 |
| `glm_1` | GLM-5.3, opencode, max | with the definitions | 107.7 min | 8 | 100 |
| `codex_nodef1` | GPT-6.1 Sol | without them | 13.8 min | 8 | 100 |
| `kimi_nodef1` | Kimi K3.1 Pro | without them | 17.9 min | 8 | 100 |
| `deepseek_nodef1` | DeepSeek V4.1 Flash | without them | 63.6 min | 6 | 50 |
| `glm_nodef1` | GLM-5.3 | without them | 147.0 min | 8 | 100 |

Each run had only the task text in its folder, Python 3.12 with NumPy and SciPy, 2 CPUs, a three-hour limit and no web
tools. One run per model of each version.

## Check an answer

Write `results.json` in the form of section 5 of the task you solved, then run

    python3 scorer/score.py <folder containing results.json> <output folder> [no_definitions]
    python3 scorer/rubric.py <output folder>

(the third argument for the second version). Python 3, standard library only.

## Check this package

    python3 check_package.py

It copies this folder to a temporary directory, checks every file against `SHA256SUMS` (and that no other file is
present), rebuilds the second version's reference, scores the eight saved answers again in the version each received and
compares the checks and points with the saved ones, and runs the controls of both versions. Its last line is
`OK: all checks passed`. Python 3, standard library only. Rebuilding the phase angles themselves takes NumPy and SciPy
(`key/reference_bilayer.py`, about five minutes) or dolfinx (`key/second_route/`).

## What needs human judgment

The scorer compares numbers. Whether a run derived, computed or recalled its formulas needs a reader; the grades in
`runs/<run>/grade/` are by Grok 4.7, and its answers on the checks a run claims are the least stable part. No expert has
reviewed the benchmark or the grades.

## Credit and status

Problem: Z. Suo and J. W. Hutchinson. Answers in `runs/<run>/submission/` were written by the AI agents named above; the
grades in `runs/<run>/grade/` by Grok 4.7; the second finite-element computation by a separate AI agent, as its report
states. The problem, the key, the scorer and the records were prepared by Teng Zhang (Syracuse University) with AI
assistants, as each record states. Corrections are welcome as issues in the MechanicsHarness repositories.

Licence: the code (the `.py` and `.sh` files) is under the MIT License (`LICENSE`); the text and data are under CC BY 4.0
(`LICENSE-CONTENT`). The source paper is not included and is not covered.
