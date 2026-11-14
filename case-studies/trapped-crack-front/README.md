# Trapped crack front: a planar crack held by tough bands under tension and two shears

A MechanicsHarness case. The case page explains the mechanics, the solution and what each agent did:
https://mechanicsharness.org/case-studies/trapped-crack-front.html

The theory comes from H. Gao and J. R. Rice (1986), "Shear stress intensity factors for a planar crack with slightly
curved front", *J. Appl. Mech.* 53, 774–778, doi:10.1115/1.3171857: the first-order change of the three stress intensity
factors along a slightly curved front of a half-plane crack (their Eq. 15). The problem built on it, a front in
equilibrium with periodic tough bands under mixed loading, is our own configuration. The agents received the problem
without the paper. The paper is not included here.

This folder lets you attempt the task, check a new answer, and check every saved answer and grade against the key.
The key is published with the case: the page and `key/` state the theory that the task withholds.

## The task in one paragraph

An infinite elastic solid (ν = 0.3) holds a planar crack x < a(z) on the plane y = 0, loaded from infinity with stress
intensity factors in the ratio K₁ : K₂ : K₃ = 1 : 0.6 : 0.8 for a straight front. Part A asks, for a prescribed wavy
front a = ā + A cos(2πz/λ), the first-order cosine and sine coefficients of the three stress intensity factors and of the
energy release rate (8 numbers). Part B asks for the front in equilibrium when the toughness varies periodically in
bands, Γ = Γ₀[1 + ε g(z/λ)]: the load level, the height of the front, the front shape and the three stress intensity
factors at 16 points each, and the largest K₃ and its position (68 numbers). A number passes within 1 % of the reference,
or 0.005 if that is larger. In the key the shear weights are (2 − 3ν)/(2 − ν) and (2 + ν)/(2 − ν), the slope couplings
2/(2 − ν) and 2(1 − ν)/(2 − ν), and the slope terms cancel in the energy release rate.

## What is here

| Path | Content |
|---|---|
| `task/PUBLIC_TASK.md` | The problem exactly as given to the runs of the first round. |
| `task/REMINDER.md`, `make_workspace.py` | The one-sentence reminder of the second round and the script that builds a run's workspace with or without it (`python3 make_workspace.py WS reminder`). |
| `key/reference.py`, `key/reference.json` | The reference: closed form from Gao and Rice's Eq. 15, all 76 numbers. |
| `key/check_direct.py` | A check of the reference that does not use its closed form for Part B: the front solved in real space and put back into the operator (NumPy, SciPy). |
| `scorer/score.py`, `scorer/rubric.py`, `scorer/RUBRIC.md` | The 76 checks of the task's rule, and the points by piece of physics (seven pieces of the response at 10 points, the front at 30, counted only when the response is complete; written after the runs, approved on 13 November 2026). |
| `scorer/check_controls.py` | Controls: the key, seven wrong first-order theories solved consistently, the tolerance edge, broken files. |
| `grader/` | The reading of each run's files by a second AI model (Grok 4.7): the instructions (`GRADING.md`, `KEY_NOTES.md`), the tools that build its read-only packet, run it and accept a grade, and their tests (`check_tools.py`). |
| `runs/<run>/` | Each run: the facts of the run (model, prompt, benchmark commit, workspace hashes, task variant), `AUDIT.md` (its tool calls and any use of the network), `submission/` (every text file it saved, as delivered), `score/`, `grade/`, and `messages.md` (its visible messages and tool steps; no private reasoning). |
| `runs/RESULTS.md` | The results of both rounds and what they show. |
| `runs/glm_1/organizer_rerun/` | Two of GLM-5.3's checks of its own elastic computation, rerun by us without network. |
| `runs/deepseek_1/deleted_working_files/` | Two working scripts that DeepSeek deleted before it finished, restored from its log. |
| `runs/recall_probe_2026-11-13/` | The question put to two models after their runs (the paper named): not part of the benchmark. |
| `runs/launch_*.sh`, `runs/make_opencode_home.py`, `runs/audit_*.py`, `runs/collect_record.py` | How the runs were launched (container, time limit, web tools off) and how their records were made. |
| `knowledge_levels/MISSED_PHYSICS.md` | The physics the three runs missed, as we read it. |
| `check_package.py` | Checks this folder (below). |
| `RELEASE_EDITS.md` | Every change made to a file for this release. |
| `SHA256SUMS` | Hashes of every file. |
| `LICENSE`, `LICENSE-CONTENT` | MIT for the code; CC BY 4.0 for the text and data. |

Not included: the full event streams of the runs and of the grader, the source paper, and the organizer's private notes.

## The runs

| Run | Model | Round | Time | Checks | Points |
|---|---|---|---|---:|---:|
| `codex_1` | GPT-6.1 Sol, Codex CLI, reasoning max | first | 5.5 min | 76 | 100 |
| `kimi_1` | Kimi K3.1 Pro, opencode, max | first | 5.9 min | 12 | 30 |
| `deepseek_1` | DeepSeek V4.1 Flash, opencode, max | first | 25.8 min | 3 | 20 |
| `glm_1` | GLM-5.3, opencode, max | first | 148.8 min | 8 | 40 |
| `kimi_r1` | Kimi K3.1 Pro | second, with the reminder | 12.4 min | 6 | 20 |
| `deepseek_r1` | DeepSeek V4.1 Flash | second, with the reminder | 16.9 min | 8 | 20 |
| `glm_r1` | GLM-5.3 | second, with the reminder | 3-hour limit, no answer | 0 | 0 |

Each run had only the task text in its folder, Python 3.12 with NumPy and SciPy, 2 CPUs, a three-hour limit and no web
tools. One run per model in each round: the table shows what happened in each run, not how often it would happen.

## Check an answer

Write `results.json` in the form of section 5 of the task, then run

    python3 scorer/score.py <folder containing results.json> <output folder>
    python3 scorer/rubric.py <output folder>

The first writes `score.json` and `score.md` (each number against the reference); the second, the points by piece of
physics. Python 3, standard library only. A close number earns no point: a piece of physics is earned only when every
number that tests it passes.

## Check this package

    python3 check_package.py

It copies this folder to a temporary directory, checks every file against `SHA256SUMS` (and that no other file is
present), rebuilds the reference from `key/reference.py`, scores the seven saved answers again and compares the checks
and points with the saved ones, and runs the scorer's controls. Its last line is `OK: all checks passed`. Python 3,
standard library only.

## What needs human judgment

The scorer compares numbers. Whether a run derived its first-order response, computed it, or recalled it needs a
reader. The grades in `runs/<run>/grade/` are by Grok 4.7, which read each run's files and answered three questions
(where the theory comes from, whether the checks the run claims are in its files, where the first error is); its
answers on checks that need arithmetic are not reliable (`runs/RESULTS.md`). No expert has reviewed the benchmark or the
grades.

## Credit and status

Theory: H. Gao and J. R. Rice. Answers in `runs/<run>/submission/` were written by the AI agents named above; the grades
in `runs/<run>/grade/` by Grok 4.7. The problem, the key, the scorer and the records were prepared by Teng Zhang
(Syracuse University) with AI assistants, as each record states. Corrections are welcome as issues in the
MechanicsHarness repositories.

Licence: the code (the `.py` and `.sh` files) is under the MIT License (`LICENSE`); the text and data are under CC BY 4.0
(`LICENSE-CONTENT`). The source paper is not included and is not covered.
