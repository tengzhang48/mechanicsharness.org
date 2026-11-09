| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: compatibility, equilibrium, nonlinear law | 8 / 8 | Uses engineering shear γ = ∇w, div τ = 0, and the correct inverse power law; the resulting p-Laplacian is correct. |
| Model: admissible radial behavior | 8 / 8 | Derives 1 + (s − 1)(1 + 1/n) = 0 from nonzero finite J, giving s = 1/(n+1) = 2/9 and the correct stress exponent. |
| Model: symmetry, traction-free faces, positive slip | 6 / 6 | Odd angular extension (2), f′(±π) = 0 as zero face traction (2), and positive A for upper-face motion in +z (2) are implemented. |
| Model: magnitude fixed by J | 8 / 8 | Uses the correct stored-work factor n/(n+1) and traction/displacement term. Reproduces I = 2.0195952773 and A = 0.0009256090617515 mm^(1−s). |
| Numerical: Δw at 0.0200 mm | 8 / 8 | Reported +0.00131813 mm; reference +0.001318126778 mm. Relative error 0.000244456%, within 2%. |
| Numerical: Δw at 0.0800 mm | 8 / 8 | Reported +0.00179369 mm; reference +0.001793693738 mm. Relative error 0.000208405%, within 2%. |
| Numerical: signed τ_rz at 0.0200 mm, 60° | 6 / 6 | Reported +38.5866 MPa; reference +38.58659686 MPa. Relative error 0.00000813155%, within 2%; sign correct. |
| Numerical: signed τ_θz at 0.0200 mm, 60° | 6 / 6 | Reported +138.011 MPa; reference +138.0106166 MPa. Relative error 0.000277781%, within 2%; sign correct. |
| Numerical: engineering-shear contour radius at 0° | 4 / 4 | Reported 0.046894 mm; reference 0.04689442305 mm. Relative error 0.000902143%, within 2%. |
| Numerical: engineering-shear contour radius at 90° | 4 / 4 | Reported 0.025066 mm; reference 0.02506612349 mm. Relative error 0.000492652%, within 2%. |
| Numerical: engineering-shear contour radius at 180° | 4 / 4 | Reported 0.013398 mm; reference 0.01339840659 mm. Relative error 0.00303459%, within 2%. |
| Completion: all seven outputs | 5 / 5 | All seven requested numbers are present, clearly identified, and reproduced by the submitted code. |
| Conventions: units | 2 / 2 | Slip and radii use mm; stresses use MPa; J uses N/mm. |
| Conventions: engineering versus tensor shear | 2 / 2 | Correctly uses γ = ∇w without a factor-of-two error, including the magnitude defining the contour. |
| Conventions: face limit | 1 / 1 | Explicitly identifies the upper-face value as the limit from inside and evaluates the angular field at π. |
| Verification: genuine requirement | 5 / 5 | Checks the traction-free boundary independently of the IVP initial conditions; also computes contour flux at four radii. |
| Verification: numerical result | 5 / 5 | Unmodified rerun gives f′(π) = 1.443 × 10^−16 and J = 0.120000000 N/mm at R = 0.02, 0.05, 0.10, and 0.30 mm. |
| Verification: interpretation | 5 / 5 | Correctly interprets the tiny face residual as agreement with zero traction and the equal contour fluxes as agreement with prescribed, path-independent J. |
| Missing fracture physics | 5 / 5 | States that initiation needs an additional toughness or cohesive criterion; the stationary constitutive field alone does not determine the advance load. |

**Total: 100 / 100** (model 30/30; numerical 40/40; completion/conventions 10/10; verification 15/15; missing physics 5/5).

The copied submission ran unchanged as `scratch/solve.py`, and the reference ran successfully with all assertions passing. Outputs are retained in `scratch/submission_run.txt` and `scratch/reference_run.txt`; all seven individual tolerance checks are in `scratch/tolerance_checks.json`. Numerical errors above use the full-precision recorded reference. No reference error was found, so an alternative corrected-reference grade is unnecessary.

(a) Scientific errors

- The description “stresses tens of MPa” in `submission/ANSWER.md`, line 104, understates the stress scale: the reported τ_θz is 138.011 MPa and the stress magnitude at the requested point is about 143.303 MPa. The supplied law gives about 126.706 MPa everywhere on the |γ| = 0.0100 contour. This ancillary wording does not change any scored result; the governing model and requested calculations are correct.

(b) Unsupported claims

- The claimed Cartesian finite-difference equilibrium residual of a few × 10^−6 (`submission/ANSWER.md`, line 92) has no implementation, step size, or numerical output in the submitted script.
- The claimed ~10^−12 agreement across DOP853, Radau, and grids (`submission/ANSWER.md`, line 102) is not demonstrated: the submitted script uses only DOP853 and contains no comparison study.
- The assertion that “the same code” ran the n = 1 limiting check (`submission/ANSWER.md`, line 87) is not demonstrated by the supplied script, which fixes n = 3.5 and contains no limiting-test branch or output. The stated analytical linear-limit formulas are correct. Full verification credit is independently supported by the reproduced boundary and contour-flux checks, so these extra undocumented claims cause no rubric deduction.

(c) Presentation problems

- The angular ODE defines p = s = 1/(n+1), whereas the usual p-Laplacian parameter is p = 1 + 1/n. The local definition makes the ODE correct, but using a different symbol would avoid confusion. All requested outputs, signs, units, and the face limit are otherwise clear.

Tolerance: **all seven requested numbers are within the stated 2% relative tolerance**, with correct signs; the largest relative error is **0.00303459%**.
