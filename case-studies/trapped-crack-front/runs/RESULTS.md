# Agent runs of version 0.3, 13 Nov 2026

One run per model. Same workspace (only `PUBLIC_TASK.md`, sha256 `c659df99…`), same prompt, same Docker image, 3-hour
limit, 2 CPUs, each model's thinking at its maximum, web tools off (see "Run conditions" in `../README.md`). The three
opencode runs record commit `fee6a92`; the Codex run records `1de5f8a`, which differs only in `runs/` and the README.
Each number below is one run: it shows what happened in that run, not how often it would happen.
A second round with one sentence of reminder added to the task (13 Nov 2026, evening) is reported in its own section
below: Kimi 20 points, DeepSeek 20, GLM 0 (no answer at the time limit).

| Run | Model | Ended | Time | Points (physics rubric) | Checks passed (A; B: level, f, κ1, κ2, κ3, max) | Coefficients used (p1, p2, p3, b2, b3, m3 against the key) | S_G |
|---|---|---|---|---|---|---|---|
| codex_1 | gpt-6.1-sol, effort max, Codex CLI 0.167.2 | by itself | 5.5 min | **100** | **76/76** (8; 2, 16, 16, 16, 16, 2) | all six as the key | 0 |
| kimi_1 | Kimi K3.1 Pro on Azure, max, opencode 2.2.6 | by itself | 5.9 min | **30** | **12/76** (3; 1, 3, 2, 1, 1, 1) | p2 = p3 = 1 and b2 = 1 (wrong); p1, b3, m3 right | −0.468 |
| deepseek_1 | DeepSeek V4.1 Flash, DeepSeek API, max, opencode 2.2.6 | by itself | 25.8 min | **20** | **3/76** (2; 1, 0, 0, 0, 0, 0) | p2 = p3 = −1 and b2 = b3 = 1 (wrong); p1, m3 right | −1.137 |
| glm_1 | GLM-5.3, Z.AI coding plan, max, opencode 2.2.6 | by itself | 148.8 min | **40** | **8/76** (5; 1, 1, 0, 0, 1, 0) | p1, p2, p3, m3 right; b2 and b3 exchanged | −2.273 |

The key: p2 = (2−3ν)/(2−ν), p3 = (2+ν)/(2−ν), b2 = 2/(2−ν), b3 = 2(1−ν)/(2−ν), m3 = 1/(1−ν). S_G is the sine
coefficient of G for the cosine front of Part A. It must vanish: G is the derivative of one energy with respect to the
front position, so the response of G to the front shape is a symmetric operator and cannot contain a slope term. This
requires b3 = (1−ν) b2, whatever the source of the coefficients.

## The physics rubric (Teng, 13 Nov 2026, after the runs)

"answer close is one thing. missing a physics is another. If one agent misses a physics, they need to lose points";
"We just need to update the rubric". `scorer/rubric.py` weighs the numbers that pass (the pass rule of the task is
unchanged) by the physics they test:

- the response, 70 points: seven pieces at 10 points, each tested by one or two numbers of Part A that isolate it and
  earned only if all of them pass: tension (C1, S1); energy release rate and load level (K10_scaled); the weight of the
  shear across the front (C2) and of the shear along it (C3); the two slope couplings (S2; S3); the energy release rate
  of the wavy front with no term in the slope (CG, SG);
- the outcome, 30 points: the equilibrium front (10) and the stress intensity factors along it (20), counted only in a
  run whose response is complete, then in proportion to the numbers that pass.

A close number earns nothing. Points: codex_1 100; glm_1 40 (both slope couplings and the energy piece missing);
kimi_1 30 (both shear weights, one slope coupling and the energy piece missing); deepseek_1 20 (only tension and the
load level). Per run: `<run>/score/rubric.md`. Controls (`scorer/check_controls.py`, part 5): the key earns 100; each
of the seven wrong theories loses exactly the pieces it misses and the whole outcome (10 to 50 points); one wrong table
value in an otherwise right answer costs 0.4 points. The weights (10 a piece, 30 for the outcome) and the rule that the outcome counts only with a complete response
were proposed by us and approved by Teng on 13 Nov 2026 ("please document. I think it is ok."); the rubric in words is
`../scorer/RUBRIC.md`. The agents were told only the pass rule; if the rubric is to be stated in the task, that is a new
version of the task text (not done).

## What each run did

- **codex_1.** Wrote the six coefficients down as "the first-order in-plane crack-front weight-function relations",
  with no source and no derivation ("No web sources or external data were used"), checked the ν = 0 limit of the slope
  terms and the cancellation of the sine term in G, solved Part B in closed form (the maximum of κ3 from the
  stationarity equation) and checked the line operator by quadrature in real space. No independent elastic computation.
- **kimi_1.** Wrote down a first-order response with the mode 1 nonlocal weight for all three modes and a slope
  coupling of 1 into mode 2 (the mode 3 coupling is right). No source, no derivation. Its eight checks are residuals of
  its own equations (about 1e-15). Its note sees the sine term in G and explains it as mode coupling.
- **deepseek_1.** Attributed to "Rice (1985)" a response in which the shear modes have the opposite sign of the mode 1
  nonlocal term ("shear modes unstable") and slope couplings of 1. Its note calls this "the main unverified assumption
  of this solution". Two of its five claimed checks (a mode 1 limit, two-dimensional limits) are not in its files.
- **glm_1.** The only run that computed the response: the crack-plane operator in Fourier space, mode 1 by an exact
  factorisation, the shear modes by a least-squares solution of the finite-part integral equation for one Fourier mode
  (27 files). Its nonlocal weights equal the key's. Its slope couplings are the key's two values exchanged, and with
  them G has a sine term. Its two ways of reading the stress intensity factors from its solution (from the traction
  ahead of the front, from the slip behind it) agree on the nonlocal weights and disagree on the mode coupling; it kept
  the first and rejected the second on the ground of an "exactly solvable inclined-front anti-plane problem" that is
  not in its files.

## Is the task statement the cause? (Teng, 13 Nov 2026: "Is our problem confusing"; "Let's figure out the real reason")

- **The statement was read as we read it.** With each run's own six coefficients, the general solver of
  `scorer/check_controls.py` reproduces all 76 numbers of kimi_1 and of deepseek_1 (largest difference 2e-9). Every
  definition, sign, normalisation, sample point and the equilibrium condition were understood as intended; only the
  coefficients differ. glm_1 begins "The task is clear (I've read it fully)". No run's messages question the wording.
- **The task can be solved from its text alone:** codex_1.
- **Recall, asked directly** (outside the benchmark, no tools, the paper named; records in
  `<review-folder>/recall_probe_2026-11-13/`): Kimi K3.1 Pro states the
  Gao-Rice formulas exactly as the key (it marks the two slope terms "informed guess"); DeepSeek V4.1 Flash states a
  different, wrong set and marks several of its coefficients "sure". So Kimi has the result and did not connect the
  problem to it; DeepSeek does not have it. GLM-5.3 was asked once and gave no answer within the 15-minute limit set
  for the question (cause not determined); gpt-6.1-sol was not asked.
- **The key.** It equals the printed Eq. 15 of Gao and Rice (1986) and Eqs. 15-17 of Lazarus (2011), and it satisfies
  b3 = (1−ν) b2. A computation that does not start from the paper now exists for the nonlocal weights only: glm_1's
  code, rerun by us in a container without network (`glm_1/organizer_rerun/`), gives p2 and p3 of the key to about
  1e-5 at ν = 0.2, 0.25, 0.3, 0.35, by both of its routes. We have not checked that code line by line. The slope
  couplings have no such computation: glm_1's two routes disagree on them, and the one it kept violates the energy
  condition above.

## Why the three runs scored so low (Teng, 13 Nov 2026: "I am also very interested in knowing why other agents fail so badly")

- **The scores overstate the distance.** 68 of the 76 numbers describe the equilibrium front, which depends on all six
  coefficients at once, so one wrong coefficient spoils almost all of Part B. By coefficients the runs stand at 6
  (codex_1), 3 (kimi_1), 2 (deepseek_1) and 4 (glm_1) of 6.
- **The count also hides how different the three answers are** (Teng, 13 Nov 2026: "Is the grading too harsh? We may
  loose some details"). Pass or fail at 1% treats a near miss like a reversed answer. Within ten times the tolerance
  (10% or 0.05) lie 63 of kimi_1's 76 numbers, 34 of glm_1's and 5 of deepseek_1's; 42 of deepseek_1's numbers have
  the wrong sign (1 of kimi_1's, 6 of glm_1's). Every run's 76 numbers follow from its own six coefficients. The task's
  rule and the scores are unchanged; the agent board (`../case_page/agent-board-desktop.svg`, built by
  `../case_page/agent_board.py`) shows the near misses as pale cells, and the case page gives the profile: coefficients
  right, numbers that follow from them, numbers within the tolerance, within ten times it, of the wrong sign.
- **kimi_1 chose in its first message and never tested the choice.** "I'll use the coplanar first-order crack-front
  perturbation in Fourier space", before any calculation (`kimi_1/messages.md`); the form it coded gives every mode the
  tension weight. It ended by itself after 5.9 of 180 minutes; its checks are residuals of its own equations. Asked
  directly, the model states the published formulas (above).
- **deepseek_1 did not have the result and chose among guesses.** Its first script compares four sign patterns of a
  tension-like term of magnitude 1/2, a second scans six; none has a weight that depends on Poisson's ratio. It kept the
  first ("Rice: mode I stable, II/III unstable") and deleted its five working scripts before finishing, although the
  task asks for the code used. They are restored from its log in `deepseek_1/deleted_working_files/`.
- **glm_1 computed the response and lost a sign** in the coupling of the shear modes; its second reading of K disagreed
  on exactly that term and was set aside (above).
- **None of the three tested its response against something it did not come from.** Part A offers such a test at no
  cost: S_G must be zero.

The case page, in the format of the fracture case library, is `../case_page/Trapped_Front_Case_Study_2026-11-13.html`
(built by `../case_page/build_page.py`; a draft for reading, not published).

## The grader on real folders (Grok 4.7; packets R1 = glm_1, R2 = kimi_1, R3 = deepseek_1, R4 = codex_1)

| Run | Where the theory comes from | Claimed checks | First error | Attempts |
|---|---|---|---|---|
| codex_1 | stated without source | all supported | none | 1 |
| kimi_1 | stated without source | all supported | theory | 1 |
| deepseek_1 | stated with source | some unsupported | theory | 2 (first rejected: summary word) |
| glm_1 | own computation | some unsupported | theory | 5 (four ended with opencode's "Duplicate tool input start") |

- The coefficient tables of the four grades agree with our reading of the method notes and the code.
- **Stable:** graded a second time, kimi_1 gets the same three answers and the same table. Over the six answers Grok
  gave for deepseek_1 (two gradings), the source and the first error are the same in six of six and the table in five.
- **Not reliable: a check that needs arithmetic.** deepseek_1 states that its κ3 equals 3.316… g(s). That is true of
  its table (ratio 3.316667 at all sixteen points, computed by us); Grok marked it contradicted in five of six answers.
  The grade on record is the one answer that marked it supported.
- **The summary of question 2.** When Grok lists a contradicted check beside unsupported ones it writes "some
  unsupported"; the form check required "contradicted" and rejected the answer. The second grading of deepseek_1 ended
  without a grade for this reason (five of five). Since 14 Nov 2026 the tool derives the summary from the listed
  checks (below); read again that way, attempts 2 to 5 of that second grading give the same three answers as the
  grade on record for the source and the first error (stated with source, theory), with the summary "contradicted"
  (kept outside git, `grading/R3b.rederived/`).
- **"Stated with source" names a paper; it does not say the paper contains the formulas.** deepseek_1 names Rice
  (1985). We have not checked that paper for the formulas it attributes to it; they are not those of Gao and Rice
  (1986).
- **glm_1's first error.** Grok places it in the change from the fixed to the local frame. Our reading: that change is
  kinematics and is right; what differs from the key is the sign of the Poisson coupling between the shear modes in
  the fixed frame (magnitude ν/(2−ν), as in the key). The grade and this reading agree on which coefficients are wrong.
- Proposed to Teng after the runs: (1) let the tool derive the summary of question 2 from the listed checks;
  (2) run each solver's check scripts ourselves, without network, and give the grader their output. (1) was applied
  first to the bilayer benchmark and, on 14 Nov 2026, here (Teng: "yes to all"; `../grader/form.json`,
  `run_grader.sh`). The grader's own word is kept in the grade as `answer_as_written`. It changes no grade that was
  accepted before; it turned the rejected answer for glm_r1 into a grade (second round). (2) is not done.

## Second round: one sentence of reminder (Teng, 13 Nov 2026)

"I think we need to just remind the AI agents to carefully check the shear modes on curved crack front". One sentence
was put into the task, before its section 1, and nothing else was changed (`../task/REMINDER.md`; task sha256
`8ba44e2d…`, first round `c659df99…`):

> **Reminder.** Carefully check the shear modes on the curved crack front.

The same three models, the same conditions, launched at 21:51 UTC from commit `d41f331`. One run per model. This is a
variant of the task: the first-round results above stand, and these are reported apart.

| Run | Ended | Time | Points (first round) | Checks passed (A; B: level, f, κ1, κ2, κ3, max) | Coefficients used | S_G |
|---|---|---|---|---|---|---|
| kimi_r1 | by itself | 12.4 min | **20** (30) | **6/76** (2; 1, 1, 0, 1, 1, 0) | p2 = 2/(2−ν), p3 = 2(1−ν)/(2−ν), b2 = 3/2, b3 = (2−ν)/2 (all four wrong); p1, m3 right | +0.758 |
| deepseek_r1 | by itself | 16.9 min | **20** (20) | **8/76** (2; 1, 1, 1, 1, 1, 1) | p2 = p3 = 1 and b2 = b3 = 1 (wrong); p1, m3 right | −1.137 |
| glm_r1 | stopped by the 3-hour limit | 177.3 min of events | **0** (40) | **0/76**: no `results.json`, no `method.md` | none formed | none |

What each run did with the reminder:

- **kimi_r1 asked for the formula again.** It started three sub-agents (opencode's `subagent` tool) and asked them
  which of two candidate sets of coefficients is the canonical one; the third was told to "derive from first principles
  or confidently recall the canonical Rice 1985/Gao-Rice equations" (`kimi_r1/messages.md`). The first answer coupled
  tension to shear, which Kimi rejected by mirror symmetry; the second answered "The canonical coplanar result is
  (B)", and the third gave the same set as the "canonical coplanar shear block". None derived the response from
  elasticity. Kimi adopted that set: shear weights 2/(2−ν) and 2(1−ν)/(2−ν).
  These two expressions are the key's slope couplings, put in the place of the weights; its slope couplings, 3/2 and
  (2−ν)/2, are not in the key. In the first round its b3 was right; now it is not, which is the piece it lost. Its
  results give S_G = +0.758, and its note takes this for physics: "The skew (sine) part is a consequence of the
  mixed-mode shear term in Irwin's relation" (`kimi_r1/submission/method.md`).
- **deepseek_r1 changed an assumption.** In the first round the two shear weights were −1; now they are +1, the tension
  weight for all three modes: the kernel "is assumed to act identically on modes I, II and III. This mode universality
  follows from the scalar Neumann/Dirichlet structure of the three near-tip fields, but it is the main external input
  of the calculation" (`deepseek_r1/submission/method.md`, section 5). The slope couplings are still the plain
  rotation of the frame. It kept its scripts this time. S_G = −1.137, as in the first round.
- **glm_r1 set out to compute the response in three dimensions and did not finish.** Its plan, in its first message:
  "represent the crack by a dislocation density on the crack plane; the exact quasistatic traction–dislocation operator
  in 2D Fourier space is diagonal in the rotated basis", with one stiffness for slip along the wave vector and another
  for slip across it, a front of finite amplitude, and the limit of small amplitude. That operator is the right physics
  (`../knowledge_levels/MISSED_PHYSICS.md`). Twenty-six messages of building and repairing the solver followed
  (`glm_r1/messages.md`). The straight crack came out "accurate to 0.07%" (its message 16). The response of the wavy
  front did not converge: in its sweeps the tension coefficient C1, which is −π, ranged from −14 to +189, and
  with a regularisation from −0.5 to −1.4; the two slope coefficients were then near +6.8 and −2.6 where the key has
  +5.91 and −3.10 (`glm_r1/submission/sweep1.log`, `sweep3.json`, `sweep4.log`). In its message 25 it planned a
  first-order solver per harmonic, the kind of formulation its first-round run had used; the limit came first. It wrote no
  answer file, so every check fails.

What the round shows, one run per model:

- All three took the reminder seriously: each reconsidered the shear modes. None answered it with the physics that is
  missing. Kimi asked its own memory again, DeepSeek replaced one assumption by another, GLM started a computation that
  did not fit into the time.
- No run tested the one condition that would have exposed its response, that the energy release rate has no term in
  the slope (S_G = 0). The two runs that answered report a nonzero S_G in their own results.
- The points did not rise: 30 to 20, 20 to 20, 40 to 0. With one run per model this is not a rate. We have no
  repetitions, so we cannot say whether a difference of one piece between two runs of a model is more than the
  variation from run to run.
- GLM's first-round run finished in 148.8 minutes with an answer; this one had no answer at the limit. A run that
  computes needs to keep a submission ready.

Grades (packets R5 = kimi_r1, R6 = glm_r1, R7 = deepseek_r1; the task text in each packet is the text that run
received):

| Run | Where the theory comes from | Claimed checks | First error | Attempts |
|---|---|---|---|---|
| kimi_r1 | stated without source | all supported | theory | 1 |
| deepseek_r1 | stated with source | all supported | theory | 2 (first not accepted) |
| glm_r1 | own computation | contradicted (the grader wrote "some unsupported") | incomplete | 5: four ended with opencode's "Duplicate tool input start"; the fifth was rejected on 14 Nov for its summary word and accepted when read again with the summary derived from its list |

"All supported" means that the checks the solver names are in its folder and do what it says. They test that its
numbers follow from its own formulas, not that the formulas are right; the grader says so in its notes on deepseek_r1.
The grade of glm_r1 is in `glm_r1/grade/` (`attempts.txt` lists the five attempts). Records: `kimi_r1/`,
`deepseek_r1/`, `glm_r1/` (made by `collect_record.py`).

## Logs

No run's log, of either round, shows a use of the network other than the model, or a path outside its workspace
(`*/AUDIT.md`). The
login copy used by codex_1 was deleted after the run. Full event streams and the grading records are outside git, in
`<run-folders>/`.
