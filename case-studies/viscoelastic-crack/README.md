# Viscoelastic crack: steady crack growth in a viscoelastic strip

A MechanicsHarness case. The case page explains the mechanics, the solution and what each agent did:
https://mechanicsharness.org/case-studies/viscoelastic-crack.html

The problem comes from Hui and Long (2026), *J. Mech. Phys. Solids*, DOI 10.1016/j.jmps.2026.106662 (main
article and supplement). The agents received the problem without the paper. The paper is not included here. The
statement is recognizable as that paper, and `task/CONTRACT.md` asks solvers to record whether they recognized it
(`FIXES.md`).

This folder lets you attempt the task, check a new answer, and check every saved answer and grade against the key.

## The task in one paragraph

An infinite strip of height 2h deforms in antiplane shear, its edges held at w = ±Δ. A crack runs along the midline at
constant speed; behind the front a cohesive zone of unknown length carries a constant shear traction τc until the opening
reaches δc. The material is a standard linear solid (modulus ratio λ = G∞/G₀ = 0.1), with a supplied free energy that
fixes the bulk dissipation; the normalized critical opening is d = 0.2. At the speeds V = 0.001, 0.1 and 10 the task asks
for the zone length ε, the load a, the energy release ratio R, the bulk dissipation ratio Γ_D/Γ₀ and stress and strain
samples, and for the range of speeds where a steady state exists. In the reference, a bounded stress at the front ties
the load to the zone length, a = (2/π) arccos(e^(−ε/2)) and R = a²/d; the critical separation fixes ε at each speed; the
zone length diverges at the terminal speed V* = 13.4235715733; and the bulk dissipation equals R − 1.

## What is here

| Path | Content |
|---|---|
| `task/PROBLEM.md`, `task/CONTRACT.md`, `task/case.json` | The problem, the submission and accuracy requirements (they announce the scorer's gates), and the input file, exactly as given to every run. |
| `key/` | The key: the reference reproduction of the published solution (`REFERENCE.md`, with its map to the paper's equations), its code (`reference_solver.py`), the values (`reference_results.json`), its checks (`verify_reference.py`, `verification_results.json`), nine values of the paper's Figure 9 read to three decimals (`published_fig9.json`), the scorer (`score_submission.py`) and the two rubrics (`SCORING.md`, frozen before the runs; `SCORING_ANALYTICAL_V2.md`, added after launch). |
| `FIXES.md` | What the audit of 25 October 2026 found, what was changed in the scorer, and the saved answers scored again. Read it before using the key. |
| `audit_2026-10-25/` | That audit, as written. |
| `runs/<run>/workspace/` | Every saved file of the three released runs, as delivered: answer, derivation, code, results and checks. `runs/<run>/records/`: the launch record and, for the DeepSeek runs, the record of their reasoning tokens. |
| `grading/R01`, `R02` | Grok 4.6's reviews of the GPT run: R01 was interrupted without a verdict; R02 is the completed review and grade. |
| `grading/R03`, `R04`, `R05` | Grok 4.6's three attempts on the Azure DeepSeek run; none delivered a final grade. R03 and R04 could use tools; R05 received a fixed packet with its tools disabled. |
| `grading/review-*` | The records of both review rounds: input freezes, runtime and delivery audits. |
| `grading/<run>-auto/` | The organizer's automatic scoring and reproduction of each run, as recorded at the time. |
| `regrade_2026-10-26/` | The scorer of 26 October run on every released answer. |
| `records/` | The pilot's process log, its protocol, the grading status, the qualitative findings, and excerpts of a review of six case pages and of its response (the parts about this case, marked inside). |
| `reference_checks/` | The review's own check of ε, a, R and V* (`viscoelastic_check.py`, output `viscoelastic_check.txt`), and the zone-length curve of the case page (`eps_curve.py`, `eps_curve.json`). |
| `diagnostics_2026-10-23/` | The diagnostics of 23 October cited by the case page (the Azure run's fit bias, exact-field controls): the script and its output. |
| `check_package.py` | Checks this folder (below). |
| `RELEASE_EDITS.md` | Every change made to a file for this release. |
| `SHA256SUMS` | Hashes of every file. |
| `LICENSE`, `LICENSE-CONTENT` | MIT for the code; CC BY 4.0 for the text and data. |

Not included: session transcripts of the agents and reviewers, the source paper, and Kimi's two attempts, which ended in
output-limit and execution failures before any mathematics was saved. Five files the study keeps outside git are not
included either: a 9 MB array file of the direct DeepSeek run (under `runs/P04_deepseek_direct_1017/workspace/dev/`) and
four records of the reviewers' tool actions.

## Check an answer

Write `results.json` as `task/CONTRACT.md` describes, then run

    python3 key/score_submission.py key/reference_results.json <your results.json> report.json

It writes every check to `report.json`, prints a summary, and exits with status 0 only if every gate passes. The gates
are those `task/CONTRACT.md` announced to the solvers: 251 items for a complete answer, 193 of them values against the
reference, three energy balances and 55 checks of counts, uniqueness and admissibility. Beside them it reports 12 checks
added after the audit of 25 October (`audit_checks_2026_10_26`): the bulk dissipation in a 5 % band, the energy balance
relative to R − 1, and the answer's own a and R against Eqs. 18 and 22. They are not counted in the gates, but they catch
an answer the gates pass: a zero bulk dissipation at V = 0.001 lies inside the announced band. Example:
`python3 key/score_submission.py key/reference_results.json runs/P01_deepseek/workspace/results.json report.json` gives
the report saved as `regrade_2026-10-26/P01_deepseek.json` (8 of 251 gates and 4 of 12 added checks fail). Python 3,
standard library only.

## Check this package

    python3 check_package.py

It copies this folder to a temporary directory, checks every file against `SHA256SUMS` (and that no other file is
present), runs the scorer on the reference values (all gates and added checks pass), reruns it on the three saved
answers and compares each report with the saved one, and runs it on deliberately wrong answers: a zero bulk dissipation,
ε just outside and inside the 0.5 % gate, the terminal speed 1/λ = 10, an a inconsistent with its own ε, a negative ε,
a stress sample just outside its gate, and malformed files. Its last line is `OK: all checks passed`. One saved report,
P04's, words one explanation as the scorer did before a correction of 26 October; the check lists it apart and compares
everything else exactly. Python 3, standard library only.

## Reproduce the key

The reference code was written for the pilot's folder layout (`agent/` for the input, `private/` for the key), and its
checks import the solver from their own folder. Run them in a copy outside this folder, with NumPy and SciPy:

    mkdir -p ../ve-check/agent ../ve-check/private
    cp task/case.json ../ve-check/agent/
    cp key/*.py key/published_fig9.json ../ve-check/private/
    python3 ../ve-check/private/reference_solver.py ../ve-check/agent/case.json ../ve-check/private/reference_results.json --energy-ny 64
    python3 ../ve-check/private/verify_reference.py

On 7 November 2026 the first command rewrote `reference_results.json` identical, byte for byte, to
`key/reference_results.json`, in about a second; the second ran its 46 checks, all passed, and wrote
`verification_results.json` equal to `key/verification_results.json` apart from its time stamp. The review's check,
`python3 reference_checks/viscoelastic_check.py`, printed the numbers of `reference_checks/viscoelastic_check.txt` (below
its two header lines), and `python3 reference_checks/eps_curve.py` rewrote `eps_curve.json` byte for byte. The
diagnostics script of 23 October records the hashes of its inputs, the source paper among them, so it does not run
without the paper; its output is kept as written.

## What needs human judgment

The scorer compares numbers only ("JSON accuracy only"). Whether an answer found the physical state, rather than a state
that passes a residual or an energy balance, needs a reader: the direct DeepSeek run chose its low-speed states partly by
energy closure and then cited energy closure as a check, and the Azure run imposed the bounded front through a fitted
coefficient that biased its zone lengths. The rubric of record is `key/SCORING.md`; `key/SCORING_ANALYTICAL_V2.md` is a
later reading, added after the GPT run had delivered.

## The runs

| Run | Agent | Model | Endpoint | Reasoning | Start (UTC) | Time used |
|---|---|---|---|---|---|---:|
| P03_gpt61_sol | GPT | gpt-6.1-sol | Azure | not recorded | 2026-10-22 14:43 | 16 min 47 s |
| P01_deepseek | DeepSeek | DeepSeek-V4.1-Flash | Azure | on | 2026-10-22 14:22 | 3 h 49 min |
| P04_deepseek_direct_1017 | DeepSeek | deepseek-flash (DeepSeek-V4.1-Flash-1017) | DeepSeek API | on | 2026-10-22 15:49 | 2 h 58 min |

Each run used OpenCode in a network-isolated container, with a six-hour budget (`runs/<run>/records/RUN.json`). The
DeepSeek runs' reasoning is recorded in `runs/<run>/records/THINKING_VERIFICATION.json`: reasoning tokens in all 226
requests of P01 and in 209 of the 246 of P04. The two DeepSeek runs differ in
endpoint, timing and route, so their results do not isolate a provider or a reasoning setting.

## Credit and status

Problem and paper: Hui and Long. Answers in `runs/` were written by the AI agents named above. The problem
statement, the reference reproduction and the scorer were prepared by the coordinating agent, Codex, from the OpenAI
model family, which is also the family of run P03 (`records/PROCESS_LOG.md`); for that reason the automatic scoring of
P03 was withheld, and Grok 4.6 replayed the scorer and reviewed its derivation against the paper. The reviews are by
Grok 4.6. The case is maintained by Teng Zhang (Syracuse University); its checks and records were prepared with AI
assistants, as each record states. No expert has reviewed the key or the grades. Corrections are welcome as issues in
the MechanicsHarness repositories.

Licence: the code (the `.py` files) is under the MIT License (`LICENSE`); the text and data are under CC BY 4.0
(`LICENSE-CONTENT`). The source paper is not included and is not covered.
