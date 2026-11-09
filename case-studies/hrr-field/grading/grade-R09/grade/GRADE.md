# Grade

| rubric item | points awarded / max | reason |
| --- | --- | --- |
| Constitutive model and inner limit | 30 / 30 | 20: selects the hardening-dominated inner field and explains why it overtakes elasticity for a singular stress field at n = 4; the later strain ansatz and code correctly use strain proportional to r^(-np). 10: correct leading plastic plane-strain closure, sigma_zz = (sigma_rr + sigma_theta_theta)/2, and effective stress squared = 3D^2/4 + 3T^2. |
| Geometry, boundary conditions, branch | 20 / 20 | 10: free faces correctly placed at theta = +/-2pi/3, with f = f' = 0. 5: correct symmetry, f(0) = 1, and positive opening amplitude. 5: obtains the reference leading tensile branch and explicitly checks finite energy, 5p = 0.7808763 < 2. The rerun also gives positive f throughout the sampled interior. The uniqueness claim is stronger than the scan establishes, but the selected branch and energy check are correct. |
| Numerical exponent p | 15 / 15 | Reported p = 0.1562 has absolute error 0.000024739714 < 0.0003. The more precise 0.15618 and 0.15617526 also satisfy the tolerance. The submitted shooting code reruns to p = 0.156175260286088. |
| Numerical stress ratio R | 15 / 15 | Reported R = 0.3733 has absolute error 0.000037686688 < 0.002. The more precise 0.37326 and 0.373262 also satisfy the tolerance. The coherent Airy-stress convention gives R = [(2-p) + f''(0)]/[(2-p)(1-p)], and the rerun gives 0.373262313311593. |
| Quantitative verification | 15 / 15 | 10: rerunning verify3.py successfully reproduces adaptive collocation with 9,439 nodes, p = 0.1561752603, R = 0.37326231, and the n = 4 crack exponent 0.2000000000. It states agreement with shooting. 5: reproduced maximum BVP residual 9.99e-12 and boundary residual 5.65e-34, with an additional successful shooting tolerance study. The inaccurate ancillary accuracy claims below do not negate these successful checks. |
| Justification and scope | 3 / 5 | Explains amplitude independence, homogeneity, and free normalization, earning the one-aspect partial-credit level. It does not explain why the stipulated monotonic proportional loading permits total deformation theory, or explicitly identify global matching as what fixes the absolute amplitude. |

**Total: 98 / 100.** No discretionary deductions beyond the rubric's explicit score levels are applied to ancillary errors.

The reference's requested values and governing equation check out. Its rerun reproduces the recorded results; the submission's explicit ODE agrees with the reference ODE to 2.67e-15 over sampled states. No corrected-reference grading is needed.

Scientific errors:

- ANSWER.md gives the elastic 240-degree wedge ratio as 0.987. The analytic symmetric free-wedge solution gives R = 0.444272010788, agreeing with the reference's independent linear-limit check.
- The first model bullet writes plastic strain proportional to sigma_e^(n-1), omitting the deviatoric-stress factor. At fixed angular stress shape the plastic strain scales as sigma_e^n; its ratio to elastic strain scales as sigma_e^(n-1). The later ansatz and implemented equations use the correct scaling.
- The auxiliary verify.py prints sigma_zz using half of [(2-p) + sigma_theta_theta], omitting f''(0). That expression gives approximately 1.69984 instead of 1.06831. ANSWER.md and solve.py's report use the correct value, so the implemented plane-strain closure receives full credit.

Unsupported claims:

- A finite grid and selected root seeds do not establish a unique symmetric eigenmode. The scan also restricts p below 0.2499; finite energy alone permits p < 2/(n+1) = 0.4. Requiring decaying displacements adds a restriction beyond the problem's stated energy condition.
- The shooting study does not make q and R identical to ten decimal places: reruns span about 5.02e-10 in q and 3.22e-10 in R. These differences are nevertheless far below the requested tolerance.
- The claimed pointwise finite-difference residual of about 2e-7 is the diagnostic's median, not its maximum. Rerunning diag3.py gives median 1.863e-7 and maximum 0.8268, with the worst error at the truncated grid endpoint. This diagnostic does not establish a uniform residual bound; the successful adaptive BVP check supplies the verification credit.
- Failure of the attempted low-order cosine fits does not establish that a global basis cannot resolve the eigenfunction.

Presentation problems:

- The answer repeats the requested values and adds inaccurate ancillary comparisons to a request for a short justification and one verification.
- Writing the problem's acceptance tolerances as +/- values can be mistaken for estimated numerical uncertainty.

**Requested numbers within tolerance: yes, both p and R, including the least precise displayed values.**

Reproduction evidence: [shooting audit](scratch/grading_checks.json), [submitted BVP output](scratch/submission/verify3.out), [submitted finite-difference diagnostic](scratch/submission/diag3.out), and [reference rerun](scratch/reference/reference_results.json).
