# Rubric — 100 points

Any sound analytical or numerical method is acceptable. Source identification,
paper terminology, and reuse of the reference code earn no points.
For numerical grading use the unrounded targets:
$\Sigma=133.730286053811\ \mathrm{MPa}$,
$L=4.08560913824553\ \mathrm{mm}$ per end.

| Criterion | Points | Full credit | Explicit partial credit |
|---|---:|---|---|
| Finite geometry and coupling | 12 | Treats the original half-slit $a$ and the two zones of length $L$ as a coupled finite crack; does not assume $L\ll a$. Equivalent continuum discretizations qualify. | 6 for correct zone geometry with an incomplete coupling formulation; 0 for replacing it entirely by a traction-free crack and plateau toughness. |
| Intrinsic front work | 7 | Accounts for $\Gamma_0$ at each advancing tip separately from connector work; e.g. $K_{\rm tip}^2/E'=\Gamma_0$. | 3 for recognizing both energy contributions but using only their total in a global crack formula; 0 if $\Gamma_0$ is discarded or doubled. |
| Elastic and local constitutive inputs | 6 | Correct plane-strain modulus (3) and constant closing traction acting against full relative opening (3). | Award the two 3-point parts separately; 1 per part for a correct stated idea with a consequential implementation error. |
| Remote stress | 25 | Within ±0.3% of target, in MPa. | 20 for error >0.3% to 1%; 12 for >1% to 5%; 5 for >5% to 15%; 0 above 15% or missing. |
| Zone length per end | 25 | Within ±0.3% of target, in mm per end. | Same 20/12/5/0 bands as stress. |
| Verification | 15 | Performs and interprets at least one meaningful check: an independent elasticity calculation, a known zero-front-work/large-crack limit, or a numerical refinement combined with physical residual checks. | 8 for checking only the imposed tip/opening residuals; 3 for naming a relevant check without executing or interpreting it; 0 for no check or a contradictory check. |
| Conventions and event | 10 | Correct full gap and per-end length/area conventions (5); identifies the root cutoff on a physically consistent monotonic branch (5). | Award each 5-point part separately; 2 for a partly correct but unclear convention or event argument. |

The numerical rows are graded independently of the method rows: unsupported
correct numbers earn their numerical points, but no unshown reasoning or checks
are inferred. If an otherwise sound method has a single arithmetic or units
error, preserve its conceptual credit and apply the numerical bands to the
reported quantities. A correctly converted SI answer is acceptable.
