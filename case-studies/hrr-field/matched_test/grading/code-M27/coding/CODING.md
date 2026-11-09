# Boundary-check coding

- **B1 existence — present_right; before.** The finite, nonzero J condition selects the admissible exponent, and the eigenvalue scan confirms that exponent for the angular field. For n = 3.5 this is consistent with the reference's singular, locally integrable-energy field. Exponent selection starts in the initial derivation, before the first ODE test (submission/PROCESS.md:2).
  - submission/ANSWER.md:49: “A finite, non-zero J therefore requires `s = 1/(n+1) = 0.222222`”
  - submission/PROCESS.md:3: “eigenvalue scan recovers s=1/(n+1)”

- **B2 limiting cases — present_right; after.** The same numerical method recovers the known linear mode-III exponent, flux normalization, and energy-release relation. The log records this as a later verification, after the main numerical calculations.
  - submission/ANSWER.md:87-88: “For `n = 1` the same code gives `s = 1/2`, `I = π`”
  - submission/PROCESS.md:6: “n=1 limit reproduces I=pi and G=K^2/(2mu)=0.120 exactly”

- **B3 applicability — present_right; after.** The answer restricts its small-strain continuum assumptions to the requested range, states the pure-power-law constitutive scope, and identifies the missing fracture criterion needed for advance. These are physically consistent with the stated stationary problem. This discussion is explicit in the final answer, whose writing follows the calculations (submission/PROCESS.md:7).
  - submission/ANSWER.md:76-77: “small strain and continuum are taken as valid in the requested range”
  - submission/ANSWER.md:110-111: “advance requires comparing it with a separately measured critical value (toughness / cohesive criterion), which these results do not contain”

- **B4 degenerate cases — absent; unclear.** Neither file examines coinciding exponents, a singularity threshold, or another parameter-dependent change in field structure. There is no supporting quote or timing to assign.

**Notes:** Timing follows the step log: finite-J exponent selection begins in the initial derivation (line 2), the n = 1 verification follows the numerical calculations (line 6), and applicability is explicit in the final answer written at line 7. Neither the n = 1 comparison nor the statement that the specified pure power law has no elastic–plastic transition checks a parameter degeneracy or structural threshold.
