# Knowledge entries: governing equations and constraints, not a recipe

Design note, 2026-10-06 (log P237). It records a discussion between Teng and Opus after the
transfer-probe pilot, together with GPT's review of that pilot. **Status: direction agreed in
discussion; the format has not been changed yet.**

## Naming (Teng, 2026-10-06, P238)

| Term | Meaning | Example |
|---|---|---|
| **framework** | Reserved for the big picture: the new way of doing mechanics research with agents | — |
| **topic** | The cross-paper level: shared governing equations, constraint checks and scope | crack-tip asymptotics (Williams, HRR, Knowles–Sternberg, notches) |
| **knowledge entry** | One paper inside a topic: its configuration, conventions, reference data, optional reference code | the Hutchinson 1968 entry |
| unit (historical) | What the three pieces tested in the probe pilot were called; kept for those records only | `corpus/knowledge_v1/P015_hutchinson_1968` |

Teng preferred "knowledge card" in spirit but not the word "card"; "knowledge entry" and "topic"
were chosen from the offered options. Below, "unit" refers only to the tested v1 pieces.

Sources: `corpus/probes/PROBE_PILOT_SUMMARY.md` (the pilot);
`corpus/reviews/2026-10-06-probe-pilot-review-gpt.md` (GPT's review, verbatim);
`corpus/batch1/HRR_CASE.md` (why the HRR unit was hard).

## Teng's points

1. "That also depends how we write the unit. We do not want the unit as a constraint."
2. "We should try to find knowledge beyond the specific problems. The framework, the
   constraints the codes should satisfy are more important. For many powerful AI agents,
   giving a specific codes may be a bad idea."

## What the pilot says (one run per cell; descriptive only)

GPT split the 18 with/without pairs by task type; Opus recomputed the split from
`corpus/probes/scores_*.json`:

| Task type | With unit | Without | Reading |
|---|---|---|---|
| T1: the paper's method, new parameters | 98.7 | 61.2 | Specific code and formulas help most here |
| T2: a variant the paper does not solve | 81.3 | 93.8 | Mixed, and the mean was lower with knowledge: 3 pairs favour it, 2 tie, 1 failed by timeout (GLM on the HRR notch). Excluding that pair: 97.6 vs 93.0 (supplementary; corrected at P239) |
| T3: applicability (valid / invalid / borderline) | 83.3 | 73.8 | Helps in 5 of 6 pairs; serious errors remain |

- **The unit as a constraint (one case).** The HRR unit's public solver has the crack angle
  built in (`solve_angular(n, plane_strain, s, normalize, n_theta)` has no wedge-angle
  parameter). On the 120° notch, GLM tried to adapt the unit's crack solver, went into a
  private function (`hrr._make_rhs`), and timed out without an answer. Without the unit it
  wrote a general solver from the governing equations and scored 98.
- **The unit as a starting point (same probe).** DeepSeek used the unit as a known limit,
  starting from the crack solution and stepping to the notch, and scored 88. Without the unit it
  used wrong radial exponents in its equilibrium equations (76).
- **Knowledge transferred without code (T1).** DeepSeek wrote its own HRR solver but took J
  from the unit's documented equation: 100 with the unit, 55 without, where its J left out the
  shear-traction work. What transferred was a verified relation and its convention, not the
  implementation.
- **Reading is not using (T3).** In P059, DeepSeek read the film unit about 50 times and still
  concluded that buckling was necessary for indentation-driven delamination.

## The direction: governing equations and constraints first, code optional

A topic and its entries should carry what is true beyond the paper's configuration, in a form
any solver's code can be checked against. The paper's specific solution becomes an instance and a source of
reference numbers.

| Layer | Contents | HRR example |
|---|---|---|
| **1. Governing equations and assumptions** | Governing equations; each assumption and where it enters the derivation; conventions; the conservation laws used and when they hold | Deformation-theory power law; Airy function; compatibility; the leading exponent follows from J only for a crack with parallel flanks, so for a notch it is an eigenvalue; normalization max σ̃ₑ = 1 and why I, K and the field must share it |
| **2. Constraints any correct solution must satisfy** (executable checks, solver-agnostic) | Residuals of the governing equations; boundary conditions; admissibility (finite energy, sign of opening); invariants and path independence where they hold; limit cases | Equilibrium and compatibility residuals on the returned field; traction-free faces; W ∝ 1/r (crack only); J consistent across several contours and between Cartesian and polar components (a necessary consistency check; agreement on a few contours does not prove path independence); n = 1 gives the Williams field; a notch with zero angle gives the crack |
| **3. Leaving the scope** | For each assumption: what changes if it is dropped, and what still holds | Notch: boundary conditions move to ±α and the exponent is an eigenvalue. Flow theory or non-proportional loading: the energy structure is lost. Finite strain: see Knowles & Sternberg (P086) |
| **4. Reference data** | The paper's numbers with their conventions, as check values | Table 1 I and E_n with normalization and constraint stated; known printed misprints |
| **5. Reference code (optional, separate)** | An implementation that passes layer 2, offered as a check, not as the route | Today's `hrr.py`, behind a note: "a crack-only reference; for other geometries write your own solver and check it with layer 2" |

**What a check can and cannot show (GPT's review, P239).** Each check states the conditions
under which it applies and its tolerance. Passing shows the implementation is consistent with
the stated model; it does not show that the model captures the mechanism of the new problem.
So entries keep two questions apart: is the solution correct for the model (checks), and is the
model adequate for the problem (scope, layer 3, and the agent's own justification). Every check
is challenged with deliberately wrong outputs before it is trusted.

Layer 2 is the new part. It turns "the constraints the codes should satisfy" into checks that
take a solver's output (arrays, fields or functions) and report pass or fail, whatever code
produced them. That gives a powerful agent freedom of method and still holds it to the
mechanics. It extends `tools/check_entry.py`, which checks an entry's prose against its own code,
to checking any solver's results against the topic's physics.

**Beyond single papers.** Layers 1–3 are naturally shared across papers. HRR (P015, P016),
Williams, Knowles & Sternberg (P086) and the notch variant are instances of one topic:
crack-tip asymptotics as an eigenvalue problem. A topic, with the papers as its entries, is where "knowledge beyond the specific problems" would live.

## How to test it (small, as GPT recommends)

- **Conditions** on the same paper, chosen to match the question (GPT's review):
  (B0) the new entry, layers 1–4 with the checks, no solver code; (B1) the same entry with the
  reference code added, so that only solver access differs (this tests Teng's point about
  code); (C) nothing. Optionally (A) today's v1 unit, which measures the whole redesign.
  Questions include cases where the model is useful and cases that need additional physics.
- **Probes:** fresh variant and applicability questions (T2 and T3 types), where code-first is
  expected to constrain and knowledge of the governing equations and checks to help; fixed budgets; 2–3 repeats per cell,
  so that differences can be told from run-to-run variation.
- **Solvers:** at least one weaker and one stronger agent, to test Teng's point that specific
  code may be a bad idea for powerful agents. The repository's rule stands: no solving with the
  family that authored the keys (Codex/GPT here).
- **Record** model choice and its justification, numerical correctness, completion,
  verification, and recognition of missing physics; an agent may change the model's
  assumptions if it says why. Also record whether the solver used the code, the governing equations or the checks (from imports and
  session logs), as was done in this pilot.

## Wording changes adopted from GPT's review

- "A verified unit usually helps" becomes "the unit helped in 12 of the 18 tested pairs",
  in `HOW_TO_HELP_AGENTS.md` and the probe summary.
- The by-task split is added to the probe summary, so that the overall mean does not hide the
  variant result.
