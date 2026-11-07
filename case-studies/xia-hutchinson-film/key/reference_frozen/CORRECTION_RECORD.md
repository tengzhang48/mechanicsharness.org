# CORRECTION_RECORD.md — post-review source-transcription correction

- Date/time: 2026-10-15T12:49–12:52Z (post-review task, outside and separate from the original
  preparation budget; coordinator has preserved the complete initial submission for review).
- Scope: single flagged item — the substrate-calibration sentence near Eq. (12) in REFERENCE.md
  Sec. 1. PROBLEM.md untouched; earlier deliverables not hidden or rewritten beyond the flagged
  sentence and this record.

## Exact before/after statement (REFERENCE.md, Section 1)

BEFORE (initial submission, REFERENCE.md lines 47–49):

> "ell is the interaction (decay) length of the crack-induced wake field; it is the length XH2000
> call l (Eq. 12), where the substrate spring is calibrated so that ell = sqrt(pi/2)·g(alpha,beta)·h
> (≈ 1.98 h with no elastic mismatch). All interaction effects enter only through H/ell."

AFTER (corrected):

> "ell is the interaction (decay) length of the crack-induced wake field; it is the length XH2000
> call l (Eq. 12), where the substrate spring is calibrated so that ell = (pi/2)·g(alpha,beta)·h
> (≈ 1.98 h with no elastic mismatch, since (pi/2)·g(0,0) = 1.976 from Eq. (1) with Beuth's
> zero-mismatch value; see CORRECTION_RECORD.md — an earlier draft mistranscribed the prefactor as
> sqrt(pi/2), an OCR radical-scope error; the PDF glyph geometry places the radical only over
> Ebar·h/k on the left side). All interaction effects enter only through H/ell."

Net change: calibration prefactor sqrt(pi/2) → (pi/2). Everything else in the sentence is retained.

## Source inspected (primary PDF, not OCR)

paper.pdf = Xia & Hutchinson, JMPS 48 (2000) 1107–1131.

- Eq. (12): file page 7 = journal page 1113, second display equation. Primary-PDF glyph geometry
  (PyMuPDF raw character/drawing extraction; the rendered clip could not be visually read by this
  model, so glyphs+rules were decoded from coordinates):
  - Radical construction (AdvMT_EX U+0081 pieces + AdvMT_RM U+0016 hook + EX top bar) spans
    x ≈ 153–176, covering only the left side; inside it: italic "E" (with macron) and "h" as
    numerator (y = 148.1), a drawn fraction rule (162.7–174.3), "k" as denominator (y = 161.7)
    → sqrt(Ebar·h/k).
  - After "=" (AdvMT_SY U+0088 at x = 178.0): glyph "p" in AdvPi2 (π) at y = 148.1 above a drawn
    rule (189.5–195.0) with "2" below at y = 161.7 — a stacked fraction pi/2 with NO radical
    spanning it — followed by italic g, "(", "a", ",", "b", ")", italic h at the common baseline
    y = 154.9. Equation number "(12)" at x ≈ 456.
  - Printed form: l = sqrt(Ebar·h/k) = (pi/2)·g(alpha,beta)·h.
- Eq. (11) (cross-check): file page 6 = journal page 1112, bottom display equation. Identical
  structure: radical only over the left side sqrt(Ebar/(kh)); right side stacked pi/2 (π above rule
  172.2–177.7) times g(a,b); number "(11)". Squaring Eq. (11) reproduces Eq. (12) exactly with the
  (pi/2) prefactor.
- Internal consistency check within the source: the sentence following Eq. (12) states l = 1.98h at
  alpha = beta = 0. With Eq. (1), G = (pi/2)·(h sigma0^2/Ebar)·g, and the zero-mismatch channeling
  value 1.976·sigma0^2·h/Ebar (Beuth 1992; Hutchinson & Suo 1992), g(0,0) = 2·1.976/pi, so
  (pi/2)·g(0,0) = 1.976 ≈ 1.98 ✓. The mistranscribed sqrt(pi/2)·g would give 1.576h, contradicting
  the paper's own 1.98h. The paper.txt OCR line ("p / 2" stacked under the radical of the left side)
  is what produced the earlier misreading.

## What changed / did not change

- Changed: the transcribed calibration prefactor in one REFERENCE.md sentence (above). No other
  occurrence of sqrt(pi/2) exists in any deliverable (grep over REFERENCE.md, STATUS.md,
  reference.json, reproduce.py, REPRODUCE.md).
- Dependent source claims reviewed: REFERENCE.md Sec. 7 ("k is calibrated (XH2000 Eq. 11) to
  reproduce the isolated-crack energy release") and Sec. 5 (isolated-crack limits citing Eqs. 10–12,
  27) remain correct as stated — they make no claim about the numeric prefactor. Nothing else found
  wrong in this narrow review.
- Requested target values: **none changed.** gamma_A, gamma_B, kI_A, kI_B at ell/H = 0.5, 1, 1.5, 2
  (nu = 0.3) depend only on H/ell and nu; the model decay length ell = sqrt(Ebar·h/k) (left side of
  Eq. 12, independently derived in REFERENCE.md Sec. 1) is unaffected; the calibration constant only
  fixes ell/h for a given elastic mismatch and never enters the normalized quantities or reproduce.py.
  No rerun needed (code and reference.json contain no occurrence of the constant; confirmed by grep).

## Status

Correction applied and recorded; task stops here as instructed. No desired formula or grading
outcome was supplied or used; the mathematical target results are preserved unless substantively
corrected — none were.
