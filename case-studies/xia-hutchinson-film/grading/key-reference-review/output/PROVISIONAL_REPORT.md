# PROVISIONAL_REPORT.md

- UTC: 2026-10-15 12:40:56 UTC
- Phase: checkpoint_report
- Role: independent scientific checker (not a replacement answer-key author)

## Derived independently (from PROBLEM.md; source used as check, not as sole computation)

Wake reduction: u1=0, u2=u(x2); eps22=u', Delta_sigma22=Ebar u', equilibrium Ebar u''=(k/h)u so ell=sqrt(Ebar h/k).

Total-traction face: sigma22=0 => Ebar u'(face)=-sigma0. Midline: u2=0, sigma12=0.

Exact strip: v=sinh(a-eta)/cosh(a), a=S/(2 ell), v'(0)=-1, v(a)=0.

Case A translating-tip: one tip per period H, far-ahead u=0, far-behind spacing-H wake.
Using PROBLEM.md Delta_Pi (authoritative, not a guessed reading of source Eq. 6):
  G_A h = -int_0^H Delta_Pi_wake dx2  =>  gamma_A = tanh(H/(2 ell)).

Case B: period 2H, one moving new tip; ahead S=2H, behind S=H.
  G_B h = (Pi_bar(2H)-Pi_bar(H))*2H  =>  gamma_B = 2 tanh(H/(2 ell))-tanh(H/ell) = 2 G'(H)-G'(2H).

Modes: reflection symmetry => KII=0; G=KI^2/E => kI=sqrt((1-nu^2) gamma).

Limits (fixed h,E,nu,k): H/ell->0, gamma_A~H/(2ell), gamma_B~(H/ell)^3/4; H/ell->inf, both gammas->1.

Source anchors read in paper.txt: Eqs. 1, 8–12 (calibration k so isolated G matches Beuth; l=sqrt(Ebar h/k)=(pi/2)g h), Eq. 6 printed DW form (ambiguous sab/eab typesetting), Eq. 27 isolated kI, Sec. 4.2 Eqs. 28/30 and identity G=2G'(H)-G'(2H). OCR garbles prefactors; identity is clear. Cosmetic/source-transcription issue, not a target-quantity error in the reference.

Analytic samples at nu=0.3 (independent of author's script; same closed forms):
ell/H=0.5: gamma_A=tanh(1), gamma_B=2tanh(1)-tanh(2)
ell/H=1: tanh(0.5), 2tanh(0.5)-tanh(1)
etc. These match the tabulated reference.json values at face-value reading of the closed forms. Independent numerical confirmation of those tables is NOT yet executed in this checker session.

## Checks actually executed this session

- Read PROBLEM.md, TASK.md, REFERENCE.md, reference.json, reproduce.py, REPRODUCE.md, STATUS.md, reproduce_output.txt, paper.txt (through Sec. 4.2).
- Did NOT replay reproduce.py yet.
- Did NOT run a separate variational/full-width strip discretization yet.
- Did NOT independently recompute collocation/FD residuals this session.
- Author's saved reproduce_output.txt claims exit-0 run (colloc max err 8.549e-15; energy residual 5.551e-16; prod mismatch 2.144e-10; C1/C2/C3 fail as expected). These are author's recorded outputs, not re-executed here. Note path in that file: wrote `.../reference/reference.json` (different folder name than current `reference_submission/`).

## Incomplete work

- Replay unchanged submitted code (copy-run only; do not edit reference_submission).
- Genuinely separate numerical deciding-quantity check (full-width strip functional + natural face tractions; not calling author's analytic functions).
- Audit live convergence and wrong-control outputs vs saved text.
- ANSWER.md, results.json, checker REPRODUCE.md, independent checking code.
- Confirm G vs E vs Ebar mode relation numerically.

## Remaining uncertainty

- Mathematical consistency of the 1D-wake energy argument looks internally consistent with PROBLEM.md; this is not physical validation of a real film (Winkler k, plane-stress sheet).
- Author numerical routes still share the 1D strip reduction (u1=0); they are not a 2D tip-field solve. Checker still owes a variational discretization that does not import the closed-form v.
- Eq. 6 typesetting vs PROBLEM.md functional: treated as notation, not erratum; problem functional is authoritative.
- No desired verdict from coordinator; status not yet assigned.
