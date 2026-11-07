# Fixes after the audit of 25 October 2026

Audit: `audit_2026-10-25/AUDIT.md` (verdict READY AFTER FIXES). The statement in `task/` is kept exactly as the agents
received it; only the key and the checker were changed, and the saved runs were checked again. The files as run are in
commit `777e5ed`.

## Statement: kept, flaws named

- (c) "On what time scale, in units of t₀, does the traction relax?" defines no measure. The e-folding time of the late
  decay (0.0852), the time at which the traction at z = h has halved (0.0953) and the time to 95 % relaxation there (0.291)
  are all readings; MP-03's (c) was graded contested.
- (a) "Describe how σ_gb(z, t) evolves" does not say what the description must contain. Some graders scored the late-time
  amplitude c₁, which is not asked for (MP-02, MP-05, MP-09).
- (b) K∞ is the edge-crack value, which can be recalled; the header allows the agent's own knowledge.

## Key (`key/KEY.md`)

- (a): states that c₁ is not asked for, so an answer without it is complete, and gives c₁ = 1.5310 for reference (audit:
  collocation up to N = 480 with Richardson extrapolation; three graders give 1.531049).
- (b): K∞ is now supported by the model's own traction-free limit, 1.121522 (audit, same method), instead of the
  handbook value alone. A recalled 1.1215 passes, so this item cannot show that the model was solved.
- (c): records which readings the 0.05–0.15 band admits and which it fails.

## Checker (`key/checks.py`)

- `lambda1_ok` added (λ₁ within 1 %). The round's count of key checks (`summary_data.py`) reads only `*_ok` flags, so it
  left λ₁ out: MP-40's 3 of 4 hid an approximate λ₁, and MP-44's single pass was a typed-in K. `lambda1_class` is kept.
- `mode1_shape_ok` band 2 % → 0.5 %: 2 % admitted the paper's Appendix B estimate 0.811. Computed runs give 0.80196.

Checked: the key passes; the Appendix B answers (λ₁ = 13.36, shape 0.811) and K = 1 fail.

## Saved runs checked again

`regrade_2026-10-26/SUMMARY.md`: `lambda1_ok` is new: true for MP-02, MP-03, MP-05, MP-09 and MP-47, false for MP-18,
MP-20, MP-40 and MP-44. No other flag changes (MP-15 and MP-32 saved no results; MP-07 has no answer). With λ₁ counted,
MP-40 has 3 of 5 and MP-44 1 of 5. The case page quotes the round's counts as they were (of 4) and names this defect.
The blind grades are as written at the time.

## Not applied

- The statement edits the audit proposes ((c) as t_½ at z = h, c₁ asked for in (a)) and the checks that go with them
  (`c1_ok`, t_½ within 1 %): they change what is asked, so they belong to a new version with new runs.
