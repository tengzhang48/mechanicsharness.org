# Boundary-check coding

- **B1 existence — present_right; after.** The angular ODE and vanishing face derivative verify that the leading field is admissible; the step log records this after running `solve.py`.

  Evidence: submission/ANSWER.md:93: “\(F'(\pi)\approx-3\times10^{-17}\), confirming the closed form and eigenvalue.”; submission/PROCESS.md:4: “face traction 0; ODE shape matches closed form.”

- **B2 limiting cases — present_right; unclear.** Substituting `n = 1` into the same closed-form angular solution gives the known linear Mode-III shape. The step log does not date this check.

  Evidence: submission/ANSWER.md:50–51: “(For \(n=1\) it reduces to the classic linear Mode-III field \(F\propto\sin(\theta/2)\).)”

- **B3 applicability — present_right; unclear.** The answer states the constitutive scope, local-field scope, and small-strain, stationary, monotonic assumptions, and excludes an additional elastic term. The step log does not establish when these scope checks occurred.

  Evidence: submission/ANSWER.md:72: “**No elastic term**: used the specified pure power law exactly”; submission/ANSWER.md:100–101: “The model returns the *local* field for a *given* \(J\); it does not know the far-field load, and the small-strain, stationary, monotonic assumptions are taken as given.”

- **B4 degenerate cases — absent; unclear.** Neither file examines coinciding exponents or a parameter threshold where the singularity appears, disappears, or changes structure. No supporting quote.

B1's explicit face-condition and ODE verification follows the solve.py run; B2 and B3 have no check-specific chronology in the step log. The n = 1 limit and the eigenvalue scan at fixed n do not examine coinciding exponents or thresholds where the singularity changes structure.
