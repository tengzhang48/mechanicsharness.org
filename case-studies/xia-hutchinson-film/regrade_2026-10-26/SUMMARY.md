# Automatic checks with the new checker (26 Oct 2026)

There was no checker; the record's numeric comparison used 1e-6 absolute (`key/KEY_FREEZE.json`; counts in
`records/NUMERICAL_EVIDENCE.json`, which covers S01, S02 and S04). `key/checks.py` counts values within 0.5 % relative and reports the 1e-6 count beside it.

| Run | Within 0.5 % | Within 1e-6 | Note |
|---|---|---|---|
| S01 (GPT-6.1) | 16 | 16 | |
| S02 (DeepSeek) | 7 | 0 | 2D FEM values, 0.2–1.6 % off |
| S03 (Kimi) | — | — | nothing delivered |
| S04 (Kimi, subscription) | 16 | 16 | `numeric_rows.json` |
| D01 (DeepSeek, thinking on) | 16 | 16 | |
| K01 (Kimi continuation) | 16 | 16 | |
