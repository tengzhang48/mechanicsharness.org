# Boundary-check coding

| Check | Status | Timing | Assessment |
|---|---|---|---|
| B1 existence | present_right | after | The reported finite-energy and J-exponent checks correctly support the admissible r^−1/2 field. |
| B2 limiting cases | present_right | after | Evaluation of the full constitutive law at decreasing r checks convergence to the linear leading limit, including angular strains beyond the measured coefficient. |
| B3 applicability | present_wrong | after | Scope and missing fracture physics are discussed, but the principal-strain calculation makes the claimed inner applicability radius incorrect. |
| B4 degenerate cases | absent | unclear | No check of a parameter endpoint or threshold changing the field's structure is shown. |

Evidence (at most two short quotes per check):

- **B1:** submission/ANSWER.md:43 — “finite energy, exponent p = 1/2 fixed by J”; submission/PROCESS.md:6 — “wrote ANSWER.md (numbers, model, verifications, screens, initiation discussion)”.
- **B2:** submission/ANSWER.md:42 — “−0.004998 (10⁻⁶), −0.005000 (10⁻⁸) → the measured −0.005000”; same line — “Same check at θ = 60° reproduces C_rr, C_θθ, C_γrθ to 6 digits at r = 10⁻⁸.”
- **B3:** submission/ANSWER.md:34 — “principal strains are ε_xx = ε_yy = a·p = 3.0857e−4·p (largest in magnitude; |ε_zz| = 2.019e−4·p)”; submission/ANSWER.md:50 — “no toughness (J_IC), critical CTOD, or critical-strain/damage criterion is supplied, and damage is explicitly excluded.”
- **B4:** No evidence.

Present checks are treated as after the main run at submission/PROCESS.md:4; B2 is the reported full-law r→0 check, not a T=E benchmark. B3's applicability screen incorrectly uses a·p instead of (a+2b)·p: ε_zz is the largest principal strain and sets the correct lower bound 0.0100 mm; no structural-degeneracy check is shown.
