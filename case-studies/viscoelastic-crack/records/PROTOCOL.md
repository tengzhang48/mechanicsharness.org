# Closed book pilot protocol

The user requested starting the benchmark on 2026-10-22. Two calibration pilots began at 14:22 UTC, with the mathematical task, key and rubric frozen beforehand in `runs/PILOT_FREEZE.json` and `runs/frozen_pilot_v1/`. Independent scientific reference review remains pending; these runs calibrate difficulty and produce provisional numeric checks rather than a final graded study.

The actual solver runtime is a pinned Docker image with Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1 and SymPy 1.14.0. Both runs use that same image, a high reasoning request, a 262144-token client context limit and 65536-token output request. DeepSeek V4.1 Flash uses OpenCode 2.0.32; Kimi K3 uses Kimi CLI 2.2.0 through the Azure FW-Kimi-K3 endpoint. The different harnesses are a confound and the pilot does not estimate a model-only difference.

The containers have no network route. A Unix socket relay allows only chat-completions requests for the assigned model. Host credentials are injected by the host relay and never mounted in a container. Read-only root filesystems, dropped capabilities and separate workspaces prevent access to the paper, key, repository and other runs. Native CLI state is ephemeral. Stored coordinator logs contain tool/status metadata and token counts, not private reasoning.

## Additional OpenAI-family comparison

At the user's subsequent request, P03 runs `gpt-6.1-sol` through OpenCode 2.0.32 and the Azure Responses API. Its settings and exact same task are frozen in `runs/P03_FREEZE.json`. It receives a full six hours from its own start. Its derivative runtime image changes only the API adapter; the scientific runtime and CPU/memory limits match the original runs. OpenCode also selects a model-specific patch tool. These API/tool differences are recorded as harness confounds. The relay permits local function tools only, with hosted tools and model delegation unavailable.

P03 is an ungraded comparison because the reference key and scorer were prepared by the OpenAI family. The monitor withholds its numeric scoring and retains its artifacts for independent review. No human waiver or seal has been created. The P01/P02 calibration scoring process remains as described below; it does not apply to P03.

## Question and primary outcome

Can an agent derive and compute the finite-strip viscoelastic cohesive-crack solution from the fully stated mathematical model, without access to the source paper? The primary outcome is a supported solution meeting all accuracy, verification, reproduction and access gates. Report the individual gates, numeric errors, review score, delivery status, elapsed work and resource cost.

Do not infer a stable model ranking, a framework benefit or general crack-tip ability from one run. Use the pilot to calibrate difficulty and repair ambiguities. A later comparison needs matched conditions and repeated runs. Choose the final repetition count from the study's estimation goal and pilot variation, not from the smallest count yielding a favorable difference.

## Inputs and isolation

Only the three files in `agent/` enter a solver workspace. Do not expose this repository, source PDFs, private reference, provenance, scorer, previous conversation, browser state, organizer memory or other runs. Use an isolated filesystem identity/container and disable network access. A separate folder under the same unrestricted account is not adequate isolation.

Supply the same pinned environment, tools, physical task and six-hour budget to compared runs. The current reproduction environment is recorded in `private/SOURCE_MANIFEST.json`; verify that the planned solver environment has the promised libraries before launch. Use one solver agent per run and keep subagents/external models unavailable. Record model and version, reasoning setting, harness version, context limits, wall-clock limit and actual termination mechanism.

If a clarification is needed, answer from the fixed physical model, log it and give the same clarification to every compared run. Do not provide method hints or target values. Preserve unsuccessful and partial answers; classify timeouts and missing delivery separately from wrong mathematical answers.

## Reference independence and review

The source solution is by Hui and Long. Codex reproduced the equations, numerical key and evaluator. The organizer reproduction has not been independently scientifically reviewed or human-sealed. Do not present its self-consistency checks as an independent solution.

The repository's independence rule restricts grading a model family on keys that family authored. Record the organizer's OpenAI family contribution. A graded OpenAI-family run requires the applicable human waiver or an independently prepared key under the repository's policy. No waiver, seal or approval has been created here. A pilot must also account for the fact that the source paper may be recognizable from model memory; ask agents to report recognition, and disclose it separately.

Before a graded launch, the user reviews the mathematical statement and a separate reviewer checks the key, branches and rubric. The agent task must remain self-contained: closed book does not justify omitting a material law, boundary condition or accuracy requirement.

## Freeze and evaluate

1. Record hashes of the exact agent bundle, source PDFs, numerical key, evaluator and rubric before the first run. A file manifest records a snapshot; it is not a human seal.
2. Execute the solver in a fresh isolated workspace. Retain command/file-access logs and exact submitted files. Keep organizer scoring outputs outside running solver workspaces.
3. Re-run the submitted entry point on the public input in a fresh environment. Run `private/score_submission.py` against the frozen key. A numeric pass is provisional until derivation, integration evidence and access are reviewed.
4. Record all gates, failures and uncertainty. Reviewers examine alternate correct solutions on their own merits. If careful submissions disagree with the reference, investigate the problem and key before attributing the gap to the agent.

The numeric evaluator deliberately uses public test values. Generality of the derivation is a review item; finite samples alone cannot prove it. Any later transfer cases, perturbed material laws or hidden field locations must be defined and disclosed as a new study version, with their own verified references and acceptance rules.

## Direct-provider addition and execution retry

At subsequent user requests, P04 uses OpenCode and the direct DeepSeek API `deepseek-flash` alias, documented on October 17 as V4.1-Flash-1017. P05 is a separate fresh Kimi CLI attempt after P02's API/client failure and a successful health check. Both reuse the original frozen mathematical inputs and get a full six hours from their own starts. P02 remains part of the reported execution record. See `runs/P04_FREEZE.json`, `runs/P05_FREEZE.json` and `runs/P05_RUNTIME_PROVENANCE.json`. P04 is a provider/revision comparison; P05 is a retry selected after observing a failure. Neither is a model-only causal comparison or a preregistered repeated-trial estimate. Future transport errors retain safe diagnostic categories, and offline replay uses the actual workspace owner's UID/GID.
