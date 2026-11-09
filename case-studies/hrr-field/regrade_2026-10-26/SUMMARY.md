# Numbers checked against the rubric bands (26 Oct 2026)

`key/check_P015.py` applies the full-credit bands of the rubrics to the numbers in each graded ANSWER.md (copied in the
audit, `audit_2026-10-25/check_P015.py`); output in `check_P015_selftest.json`. The 100-point totals are as graded.

| Task | Run | Grade | Numbers | All at full credit |
|---|---|---|---|---|
| T1 | with-glm | R11 | p, H, C full | yes |
| T1 | with-dsaz | R35 | p, H, C full | yes |
| T1 | without-dsaz | R31 | p full; H 2.15092 and C 0.29178 fail (graded 55) | no |
| T1 | without-glm | R05 | no answer (timeout) | — |
| T2 | with-dsaz | R32 | p, R full | yes |
| T2 | with-glm | R22 | no answer (timeout) | — |
| T2 | without-glm | R09 | p, R full | yes |
| T2 | without-dsaz | R23 | p 0.08764 and R 0.3330 fail (graded 76) | no |
