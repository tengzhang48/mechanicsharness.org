# Harness notes: how each agent was run, and three corrections (log P393–P395)

> **Excerpt for the crack-bridging case, made on 8 November 2026.** Kept, unchanged: the title and the section on how each agent ran. Left out: the corrections and the other sections.

## How each agent ran in this round (`launch.sh`)

| Agent | Harness | Command | Reasoning in the run logs |
|---|---|---|---|
| DeepSeek-V4.1-Flash (Azure) | Cline 3.0.57, JSON mode | `cline -P openai-compatible -m DeepSeek-V4.1-Flash --json` | **none:** 0 reasoning chunks in all 12 runs |
| Kimi K3 (Azure, FW-Kimi-K3) | Cline 3.0.57, JSON mode | `cline -P openai-compatible -m FW-Kimi-K3 --json` | yes: e.g. 43k reasoning chunks in MP-22 |
| GLM-5.3 (zai coding plan) | opencode 1.17.18 | `opencode run --auto -m zai-coding-plan/glm-5.3 --format json` | yes: 90k–150k reasoning tokens per B run |
| GPT-6.1-sol (Azure) | opencode 1.17.18 | `opencode run --auto -m azure/gpt-6.1-sol --format json` | yes: 5k–11k in the main thread, plus its subagents |
| Codex (grader and reviewer only) | codex 0.160 | `codex exec`, gpt-6.1-sol, reasoning effort "max" | — |

- **Cline's provider settings** (`~/.cline/data/settings/providers.json`) have `"reasoning":
  {"enabled": false}` for the openai-compatible provider. No `--thinking` flag was passed.
- **K3 reasoned anyway,** by its own default. DeepSeek did not.
