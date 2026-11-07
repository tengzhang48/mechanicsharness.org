# Milestone 2: a gap in AI agents' mechanics, and a knowledge entry that closes it

> **Excerpt for the grain-boundary wedge case, made on 7 November 2026.** Kept, unchanged: the title and section 6. Left out: the other sections.

## 6. Corrections to our own process (`HARNESS_NOTES.md`)

- **DeepSeek first ran with reasoning off.** Cline needs `--thinking`, and the saved provider
  setting disables reasoning. All its first-round results were superseded by the rerun.
- **Three GLM failures were the provider's 5-hour usage limit,** not the agent. Two earlier
  statements to the contrary were corrected.
- **Cline JSON mode ends a run on a turn without a tool call, and on agentic compaction.**
  - Runs stopped this way: MP-07 (twice), MP-15, MP-20, MP-42, MP-47, MP-72.
  - Four were recovered by one restart each: MP-20, MP-42, MP-47, MP-72.
- **B's statement did not fix the order of limits.** The key now grades the outer value 1 as
  "partly".
- **The blind grader** (Codex, gpt-6.1-sol) is the same model family as the GPT solver.
