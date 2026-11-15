# Fracture Mechanics 101 — two audiences, shared physics

The two editions serve different purposes:

| Edition | Location | Intended use |
|---|---|---|
| Human, illustrated | [Fracture Mechanics 101, topics 0–19](../case-studies/fracture-101.html) and [mixed-mode and trapping, topics 20–21](../case-studies/fracture-101-curved-fronts.html) | Conceptual explanation, figures, worked examples, self-check exercises, study links |
| Agent, concise Markdown | [fracture-101.md](fracture-101.md) | Isolated knowledge input about transferable mechanics; no case-specific coefficients, worked case answers, keys or grader metadata |
| Optional procedural intervention | [verification-workflow.md](verification-workflow.md) | Separately test whether explicit verification guidance improves reliability |

**Status:** Candidate educational materials. The new editions and their impact have not been evaluated in matched agent runs.

## Rationale

A human learner benefits from diagrams, layered reasoning, exercises and revealing the worked calculation on request. An AI agent benefits from explicit scope, definitions, assumptions, constitutive regimes and applicability constraints in a versioned, machine-readable file. The scientific core should agree; the presentations need not be identical.

A public tutorial is **not** automatically a fair agent benchmark input. The current human page uses published numerical examples and links to solved case studies. An agent reading that page during a graded task could obtain case-relevant information without doing the intended mechanics. Do not use human HTML or the full website as the agent's scored knowledge source.

## Controlled experiment

Keep the original task unchanged and define separate arms:

- **A — Task only:** no knowledge file.
- **B — Task + physics:** supply exactly `fracture-101.md`.
- **C — Task + physics + verification:** supply exactly `fracture-101.md` and `verification-workflow.md`.

Run the same models with the same tool policy, compute/time budget and scored quantities; repeat independent runs for each arm. If the task is at the ceiling without help, construct a fresh *transfer variant* whose parameters and answer were not incorporated into either edition.

Before opening any solver session:

1. Freeze the benchmark task, any transfer variant, each knowledge file, evaluation criteria and a list of known acceptable alternative approaches. Record hashes of **exact delivered bytes**.
2. Check that the agent's workspace contains only the task and its assigned files; the scorer and independent reference stay outside.
3. Check for contamination: no special target coefficients, numerical keys, source-paper-derived tables, prior submissions, or references to the specific benchmark's answer should appear in the delivered knowledge. Check named examples as well as numbers.
4. Predefine how a run is counted (success/failure/partial), what to do on infrastructure failure, and how to report numerical accuracy *separately* from correct physics and independent verification.
5. Report both overall outcomes and the first material modeling mistake, including uncertainty from small run counts.

The previous experiments with paper-specific knowledge units and with task reminders are historical baselines, **not** the outcome of testing this agent edition.

## Edit policy

- Treat `fracture-101.md` as versioned benchmark input. Do not change it mid-round.
- A physics correction requires a new version and a recorded rationale; old results continue to refer to the old hash.
- Do not silently add `verification-workflow.md` to arm B. Without this separation, a score change cannot be attributed to physics knowledge versus explicit verification behavior.
- Do not import text from the case keys or grader notes into either delivered agent file. Use neutral illustrations in human pages and disclose any relevant overlap.
- No agent benchmark run is represented as completed by this documentation change alone.

The human extension is a companion page rather than a replacement of the long existing human tutorial, so links from the released case studies continue to work.
