# Harness notes: how each agent was run, and three corrections (log P393–P395)

> **Excerpt for the grain-boundary wedge case, made on 7 November 2026.** Kept, unchanged: the title and correction 2. Left out: the other sections.

## Correction 2: three GLM failures were the provider's usage limit, not the agent

The opencode log (`~/.local/share/opencode/log/opencode.log`) shows, from 15:55 to 17:03 UTC,
repeated `"AI_APICallError: Usage limit reached for 5 hour"` for zai-coding-plan/glm-5.3.

| Run | Problem | What we recorded before | What it was |
|---|---|---|---|
| MP-32 | E | a long thinking step, then a stall, then a watchdog stop | the usage limit (opencode retrying with back-off) |
| MP-35 | D | timeout with no answer, counted as an agent failure | the usage limit |
| MP-34 | D | timeout, with an answer | the usage limit during the run |

- **Wrong earlier statements, now corrected:**
  - P372: "MP-32 ... still thinking";
  - the round tally: "MP-35 timeout without answer (agent)".
- **The watchdog's "no CPU, no reads" record fits an opencode process sleeping in back-off.**
- **Not affected:** GLM's other runs (graded) finished before 15:55. The B runs MP-06, 10 and 21
  ended at 15:22, 14:56 and 15:41.
