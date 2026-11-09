| rubric item | points awarded / max | reason |
| --- | --- | --- |
| Constitutive model and inner limit | 30 / 30 | 20: correctly selects the hardening-dominated inner limit and explains that strains proportional to effective stress to the fourth power dominate elasticity. 10: correctly derives plastic plane strain, `sigma_zz = (sigma_rr + sigma_tt)/2`; `disp.py:20–29` implements the effective strain/stress inversion with tensor shear and the correct factors. The arbitrary positive constitutive scale cancels from the eigenvalue and stress ratio. |
| Geometry, boundary conditions, branch | 20 / 20 | 10: free tractions are imposed at the correct faces, theta = +/-2pi/3. 5: even/odd opening symmetry, an amplitude gauge, and positive opening stress are specified. 5: explicitly checks the correct finite-energy condition, `5p < 2`, excludes its computed infinite-energy branch, and follows the linear tensile branch by continuation. These are credited branch-selection criteria; the equilibrium error below invalidates the resulting eigenvalues. |
| Numerical exponent p | 4 / 15 | The most precise reported value, 0.08763859, differs from 0.156175260286 by 0.068536670286, exceeding both 0.0003 and the partial-credit tolerance 0.003. Award the nonlinear eigenvalue method level for a justified displacement BVP and continuation attempt. The incorrect equilibrium equations prevent an accurate result. |
| Numerical stress ratio R | 4 / 15 | The most precise reported value, 0.33301436, differs from 0.373262313312 by 0.040247953312, exceeding both 0.002 and the partial-credit tolerance 0.02. Award the expression-on-an-imperfect-solution level: the code evaluates the correct same-radius ratio `S_rr(0)/S_tt(0)` with a coherent stress convention. |
| Quantitative verification | 15 / 15 | 10: a relevant analytic limiting check is actually supported. Rerunning the linear notch gives p = 0.384268905147 versus the analytic 0.384268940509, an error of 3.54e-8, consistent with the stated agreement. 5: all four `conv.py` runs reproduce p = 0.08763859 and R = 0.33301436 at eight decimals and quantify residuals of 4.8e-9 to 6.2e-9. This earns refinement/residual credit, although every run returns status = 1 and the residuals concern the incorrect equations. The claimed successful crack continuation is not credited. |
| Justification and scope | 3 / 5 | Connects remote y-tension to symmetric tensile opening and performs the dimensionless calculation with an amplitude gauge. It does not explain the role of monotonic proportional history in the deformation-law model or explicitly recognize that global matching fixes the unrequested amplitude. |

**Total: 76 / 100.**

(a) **Scientific errors**

- The decisive error is the stress radial exponent in equilibrium. With the submission's notation, stress scales as `r^-p`, strain as `r^(-np)`, and displacement as `r^(1-np)`. The correct angular equilibrium equations are `(1-p)S_rr - S_tt + S_rt' = 0` and `(2-p)S_rt + S_tt' = 0`. The answer and `disp.py:57–58` instead use `np` in both stress equations; its pressure-scaling comment is also inconsistent with the constitutive stress scaling. At n = 1 this error disappears, so the successful linear check cannot detect it.
- Consequently, the n = 4 field and the identification of its roots as the physical notch eigenvalues are incorrect. If the stated erroneous radial equation were satisfied exactly, the correct radial equation would still have residual `(n-1)p*S_rr`. Using the answer's positive center stress, this is approximately 0.141054, rather than a small numerical residual. The correct finite-energy inequality alone does not repair equilibrium.
- Agreement across grids does not establish accuracy against the requested physical solution. Both final numbers fail even the rubric's wider numerical tolerances.

(b) **Unsupported claims**

- The fields are repeatedly called converged. Actual `conv.py` output gives status = 1 (maximum mesh nodes exceeded) for every row; requested tolerances of 1e-10 through 1e-12 were not achieved. The stable reported digits are reproducible, but successful convergence at those tolerances is not.
- The claimed continuation residual of at most 5e-9 “throughout” is contradicted by the submitted `contn2.py` rerun: its initial n = 1 residual is 7.41e-9. Subsequent continuation steps also return status = 1.
- A successful numerical n = 1 crack limit is not demonstrated by the provided angle-continuation code. Rerunning `conta.py 1.0` ends at beta = 0.998pi with p approximately 0.499080, status = 1, and residual approximately 1.09e2. The analytic crack limit is valid, but this calculation does not verify it reliably.
- The assertion that the two computed roots constitute the physical mode-I spectrum is unsupported because the underlying n = 4 equilibrium equations are wrong.

(c) **Presentation problems**

- The refinement table omits solver status and achieved residuals, obscuring the distinction between requested tolerance and actual convergence.
- Calling the “inplane field” deviatoric is ambiguous: the total stress includes the pressure term that the code correctly retains. The displayed center stresses should also be labeled normalized angular stress components.

Reference audit: the rerun in `scratch/reference/` exactly reproduces `reference_results.json`, including the independent shooting and limiting checks. No reference error was found; no alternative grading is warranted. Rerun evidence is saved in `scratch/submission/conv.log`, `contn2.log`, `initial_checks.json`, and `linear_angle_continuation.log`.

**Tolerance verdict:** neither requested number is within tolerance: absolute errors are 0.068536670286 for p (allowed 0.0003) and 0.040247953312 for R (allowed 0.002).
