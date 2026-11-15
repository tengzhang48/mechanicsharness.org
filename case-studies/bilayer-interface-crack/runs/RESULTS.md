# Agent runs of the bilayer problem (version 0.2), 14 Nov 2026

Teng, 14 Nov 2026: "please test the bilayer problem with the agents, including codex". One run per model. Same
workspace (only `PUBLIC_TASK.md`, sha256 `21c0749a…`, no figure), same prompt, same Docker image, 3-hour limit that
the agent is not told, 2 CPUs, each model's thinking at its maximum, web tools off (see "Run conditions" in
`../README.md`). All four record commit `3a37a37`. Launched between 02:40:52 and 02:41:54 UTC. Each number below is one
run: it shows what happened in that run, not how often it would happen.

**All four runs passed all eight checks.**

The variant without the interface-crack definitions (14 Nov 2026) is reported in its own section below: three runs
passed all eight checks, DeepSeek six.

| Run | Model | Ended | Time | Points | Checks | ψ, force (key 52.249°) | ψ, moment (key 4.279°) | ψ, both (key 30.433°) |
|---|---|---|---|---|---|---|---|---|
| codex_1 | gpt-6.1-sol, effort max, Codex CLI 0.167.2 | by itself | 11.6 min | **100** | **8/8** | 52.2494° (+0.000°) | 4.2787° (+0.000°) | 30.4333° (+0.000°) |
| kimi_1 | Kimi K3.1 Pro on Azure, max, opencode 2.2.6 | by itself | 8.5 min | **100** | **8/8** | 52.15° (−0.10°) | 4.24° (−0.04°) | 30.36° (−0.07°) |
| deepseek_1 | DeepSeek V4.1 Flash, DeepSeek API, max, opencode 2.2.6 | by itself | 41.6 min | **100** | **8/8** | 52.27° (+0.02°) | 4.25° (−0.03°) | 30.42° (−0.01°) |
| glm_1 | GLM-5.3, Z.AI coding plan, max, opencode 2.2.6 | by itself | 107.7 min | **100** | **8/8** | 52.18° (−0.07°) | 4.34° (+0.06°) | 30.48° (+0.05°) |

The five closed-form numbers (ε, the two energy release rates, the two magnitudes of K) are right to all digits in
every run. The check on a phase angle allows 0.5°. The points follow the rubric in `../scorer/rubric.py`, approved
on 14 Nov 2026 after this round; with every piece passed, the choice of weights does not matter for this round.

## What each run did

Every run obtained ε from the formula of the task, the energy release rates from the beam energies of both layers
far behind the tip, the magnitudes of K from G with the two-material relation, and the phase angles from a
finite-element solution of its own, written with NumPy and SciPy in the run. None took a phase angle from the
literature.

| Run | Finite elements | How K was read from the solution | Its second check of the phase |
|---|---|---|---|
| codex_1 | 9-node quadrilaterals, rings of elements collapsed onto the tip, three meshes up to 103 thousand nodes | crack opening behind the tip, with the two-material formula | tractions ahead of the tip; the exact difference of the two phase angles (below) |
| kimi_1 | 4-node quadrilaterals graded to the tip, three meshes up to 259 thousand unknowns | crack opening, fitted with further terms of the tip expansion | tractions ahead of the tip; it reported "central" values between its fits and stated an uncertainty of about 0.2° |
| deepseek_1 | elements of order 2 and 3, four meshes up to 290 thousand unknowns | crack opening | tractions ahead of the tip; an interaction integral with a crack-tip field it built itself; the case of both loads solved directly, the sum of the two K as the check |
| glm_1 | 4-node quadrilaterals graded to the tip | tractions ahead of the tip, as the task defines K | the exact difference of the two phase angles; two one-material cases with known K; the case of both loads solved directly. A second route (a Fourier solution of the layers with collocation) was built, showed a bias of 2 to 10% in the magnitude of K, and was set aside |

**The exact difference of the two phase angles.** K is linear in the loads and G quadratic, so the energy alone fixes the
size of ψ_A − ψ_B, 47.9705°; its sign follows from the computed factors. codex_1 and glm_1 state this and test their numbers against it; kimi_1 and deepseek_1 do not
mention it. The reported angles miss it by 0.0002° (codex_1), 0.06° (kimi_1), 0.05° (deepseek_1) and 0.13° (glm_1).

## Where the interface-crack formulas came from (Teng: "did the agents write from memory, or derive, or search")

No run searched: web tools were off, no log shows a use of the network, and the tools used were reading, writing and
running files (`*/AUDIT.md`). No run names a paper for a formula. The task gives ε, the definition of K by the
tractions and the definition of ψ. It does not give the relation between G and the magnitude of K, nor a formula for
the crack opening; the runs supplied these.

| Run | The formulas the task did not give | How the run tested them |
|---|---|---|
| codex_1 | both written directly; the relation between G and K is in its first calculation, before any solver existed | the magnitude of K from the crack opening equals the magnitude from the energy within 0.001%; K from the tractions ahead of the tip, which needs only the task's definition, gives the same phases within about 0.06° |
| kimi_1 | both written directly | it then set up the two-material tip problem itself (three scripts) to compare with its opening formula; the comparison did not close and it deleted the scripts (they are in its event stream); the tractions ahead of the tip give phases within about 0.15° of those from the crack opening |
| deepseek_1 | both written directly | it solved the two-material tip problem numerically (exponent 1/2 + iε from its own determinant) and compared: its opening factor against the formula to 1 part in a million, the energy of the unit field against the relation between G and K to the same accuracy (seen in its event stream; the outputs were not saved in its folder) |
| glm_1 | the relation between G and K written as "the standard interface result"; no crack-opening formula used | at ε = 0 by a one-material test case; at the task's ε through the agreement of the magnitude of K from its tractions with the energy value (0.08% and 0.56%) |

We see what each agent wrote and ran, not its private reasoning: "written directly" means that the formula appears
without a derivation and without a source before it is used. glm_1 is the only run that names the literature of the
problem, as something it did not have: "an independent reference (e.g. tabulated Suo-Hutchinson universal functions)
was not available offline".

## What the round shows

- **The problem as posed is within reach of all four agents.** The task hands over the physics that is special to
  an interface crack (ε, the definition of K with its oscillating factor, the reference length). What is left is
  energy accounting, two standard relations that every run knew, and a two-dimensional linear finite-element
  computation, which our own code does to 0.14° in one second (`../key/container_timing.log`).
- **It is the opposite of the trapped front** (`../../trapped_front_mixed/runs/RESULTS.md`). There the missing piece
  was physics that no cheap computation supplies, and three of four runs replaced it by a remembered form: 30, 20 and
  40 points for Kimi, DeepSeek and GLM. Here the missing piece only had to be computed, and the same three models
  computed it.
- **Four more computations of the phase angle, none of them ours.** They give 52.15° to 52.27° under the force; our
  two routes give 52.249° and 52.268°. The entry of the source paper's printed table at these parameters is 51.4°
  (see "What is not settled" in `../README.md`). The runs did not have the paper, and none refers to the table's
  value.
- One run per model: not a rate.

Teng, 14 Nov 2026, on the first three results: "Why this problem is simple"; then "How about we remove the
interface-crack physics". That variant, approved the same day, is `../task/VARIANT_NO_DEFINITIONS.md`; its runs
(`*_nodef1`) are reported in their own section.

## The variant without the interface-crack definitions (14 Nov 2026)

Teng, after the first round: "Why this problem is simple"; "How about we remove the interface-crack physics"; the
wording approved the same day ("yes to all"). The variant (`../task/VARIANT_NO_DEFINITIONS.md`, sha256 `74043e0d…`)
keeps the bilayer and the loading word for word and leaves out the formula for ε, the definition of K with its
oscillating factor and the reference length. It asks the same quantities in other words: the energy release rates, the
amplitude of the crack-tip traction (= |K|), the phase of the crack-tip field at r = h (= ψ), and the change of that phase
when r increases tenfold (= ε ln 10, −8.514°). Same conditions as the first round; launched 13:51:54 to 13:52:56 UTC
from commit `681a64b`. One run per model; reported apart from the first round, whose results stand.

| Run | Ended | Time | Points | Checks | change per decade (key −8.514°) | phase at h, force (key 52.249°) | phase at h, moment (key 4.279°) | both (key 30.433°) |
|---|---|---|---|---|---|---|---|---|
| codex_nodef1 | by itself | 13.8 min | **100** | **8/8** | −8.5136° | 52.25° | 4.29° | 30.44° |
| kimi_nodef1 | by itself | 17.9 min | **100** | **8/8** | −8.5136° | 52.29° | 4.34° | 30.48° |
| deepseek_nodef1 | by itself | 63.6 min | **50** | **6/8** | −8.5136° | 51.6° (−0.65°) | 5.4° (+1.12°) | 30.6° |
| glm_nodef1 | by itself | 147.0 min | **100** | **8/8** | −8.5136° | 52.27° | 4.25° | 30.43° |

**What removing the definitions changed: nothing that the scores can see.** Every run wrote the crack-tip field of an
interface crack, T = K r^(iε)/√(2πr), and the bimaterial constant in the Dundurs form ε = ln[(1 − β)/(1 + β)]/2π, without
a derivation and without a source, before using them, and every run got the change per decade exactly. Codex, Kimi
and DeepSeek wrote the relation between G and the amplitude with cosh²(πε) the same way; GLM gives it under the heading
"G–K relation (derived and verified)", from a crack-tip field it set up and checked in its own code
(`code/eigenpair.py`), and calls its energy release rates "the Suo–Hutchinson 'beam' values", the only run to name
the source. Each run then tested what it had written
against something of its own: Codex fitted the slope of the phase against ln r in its finite-element solution
"without imposing its value" (ε = −0.06451 and −0.06453), which "independently checks the oscillation sign and
size"; Kimi compared the magnitude from its crack-opening
formula with the one from the energy, and one material with a known K; DeepSeek checked that the determinant of the
two-material tip problem "vanishes at `s = -1/2 +/- i*eps`"; GLM measured ε from the slope of its solution and checked
the relation between G and the amplitude by a J integral of the crack-tip field it built (5 × 10⁻⁷). None took a phase
angle from the literature.

**The one miss is numerical.** DeepSeek's two phases are 0.65° and 1.12° off; the grader places the first error under
"numerics": by its reading of the code, the phase averaged over 3 × 10⁻⁵ < r/H < 10⁻³ (6 × 10⁻⁵ < r/h < 2 × 10⁻³; the
run's unit of length is H) rather than extrapolated, on a model reaching 3H = 6h behind and ahead of the tip.
Its two phases differ by 46.2° against the exact 47.97° that follows from the energy; Codex states that check in its
note (its difference 47.96°), the other three do not. The phase under both loads passes (30.6°) although it is built
from the two that fail.

**A search outside the folder.** DeepSeek ran one command that searched the container's system folders for files that
mention "suo", "interface crack" or "hutchinson" (`grep -ril … /usr/share /opt /usr/local/lib/python3.11`). It found
only unrelated system files (certificates, a changelog, a character table) and opened none of them. The prompt said
"Work only inside this folder and do not use the web"; no web was used (`deepseek_nodef1/AUDIT.md`, paths outside the
workspace). The other paths it named outside the workspace are the output files of its own background commands.

**What the round suggests.** The interface-crack relations are carried by all four models: taking them out of the task
does not make the problem harder for them. Before the runs I expected this (`../README.md`, the variant's section); one
run per model does not make it a rate. Whether this is recall from training data cannot be seen from the files; what
can be seen is that the formulas appear at once, in standard notation, before any computation, and are then checked.

| Run | Grade (Grok 4.7): where the phases come from | Claimed checks | First error |
|---|---|---|---|
| codex_nodef1 | own computation | some unsupported (7 supported, 1 unsupported, 1 not decidable) | none |
| kimi_nodef1 | own computation | some unsupported (7 supported, 1 unsupported) | none |
| deepseek_nodef1 | own computation | some unsupported (2 supported, 3 unsupported, 3 not decidable) | numerics |
| glm_nodef1 | own computation | contradicted (3 supported, 6 unsupported, 2 contradicted; the grader wrote "some unsupported") | none |

Each was graded at the first attempt (packets B5 to B8). glm_nodef1's two contradicted checks: its note states that
two ways of reading the amplitude agree to 0.1 %; its own `final_audit.txt` shows 0.16 %.

## The grader (Grok 4.7; packets B1 = kimi_1, B2 = codex_1, B3 = deepseek_1, B4 = glm_1)

| Run | Where the phase angles come from | Claimed checks | First error | Attempts |
|---|---|---|---|---|
| codex_1 | own computation | all supported (8) | none | 2 |
| kimi_1 | own computation | all supported (5) | none | 1, after ten failed attempts on the packet without the bundled file |
| deepseek_1 | own computation | some unsupported: 3 unsupported, 11 not decidable from the folder | none | 1 |
| glm_1 | own computation | all supported (9) | none | 1 |

- deepseek_1's note names a patch test and a one-material test case whose code is not in its folder, and for most of
  its other checks the code is there but no output was saved. Its event stream shows that several of them were run.
  On the trapped front this model deleted its working scripts; here it kept the solver and left out the outputs.
- On kimi_1's folder of twenty files the grader's first ten attempts ended with opencode's error "Duplicate tool input
  start". The packet now also holds every text file of the submission in one file, and the folder was graded at the first attempt
  (see "The grader" in `../README.md`). The failed attempts are kept outside git.
- For deepseek_1 the grader wrote the summary "cannot tell" beside a list that holds three unsupported checks; the
  tool derived "some unsupported" from the list, as set in `../grader/form.json`.

## Logs

No run's log shows a use of the network other than the model (`*/AUDIT.md`). One run named paths outside its
workspace: deepseek_nodef1, once to search the container's system folders (section on the variant) and otherwise for the
output files of its own background commands. The login copies used by codex_1 and codex_nodef1 were deleted after the
runs. Full event streams and the grading records are outside git, in
`<run-folders>/`. Per run here: `<run>/` with the facts of the run, the audit, the score,
the grade, the submission and the agent's visible messages.
