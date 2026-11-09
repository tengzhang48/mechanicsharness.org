| rubric item | points awarded / max | reason |
| --- | ---: | --- |
| Stress exponent — numerical value | 10 / 10 | Submitted p=0.1333333 differs from 0.133333333333 by 3.3333e-8, within 1e-5. |
| Stress exponent — contour-flux scaling | 5 / 5 | Correctly connects the constant nonzero J to W proportional to r^-1, giving p=1/(n+1), strain proportional to r^(p-1), and displacement proportional to r^p. |
| Mechanics — equilibrium and compatibility | 5 / 10 | The Airy stresses satisfy equilibrium. The reconciliation paragraph and redo.py correctly set up the shear compatibility relation; redo.py's RHS agrees with the reference at a diagnostic state to 8.88e-16. However, the final ANSWER.py solves a different, incompatible equation, so the stipulated angular problem is not solved. Award the explicit setup-only tier. |
| Mechanics — plane-strain closure and constitutive factors | 10 / 10 | Correct closure is shown in the code: s_zz=0, hence sigma_zz=(sigma_rr+sigma_theta_theta)/2. The final computation uses Q^2=3(A-B)^2/4+3(sigma_rtheta)^2, normal strains 3q(A-B)/4 and its negative, and tensor shear strain 3q sigma_rtheta/2. No elastic term is added. |
| Mechanics — symmetry, free faces, normalization | 3 / 5 | Credit the correct phi'(0)=0 condition, zero crack-face shear traction, and explicit max Q=1 normalization (one point each). The second symmetry condition phi'''(0)=0 and the normal free-face condition sigma_theta_theta(pi)=0 are absent; sigma_rr(pi)=0 is imposed instead. |
| H — numerical value/procedure | 0 / 20 | H=2.15092 has absolute error 0.257350274787 relative to 1.893569725213. It misses both 0.005 and 0.05 tiers. The wrong boundary conditions, incompatible final ODE, and omitted contour term prevent the sound-procedure alternative from applying. |
| H — calibration by the stated J | 2 / 5 | Correctly recognizes J=alpha epsilon_0 sigma_0 K^(n+1) I_n and H=B(0)/I_n^p, with the stipulated dimensional prefactor. But I_n=7.7687 drops required shear-traction work and is reproduced by a Cartesian check with incorrect shear parity. Award calibration recognition only: the contour normalization constant is mishandled. |
| C — numerical value/integral | 0 / 15 | C=0.29178 has absolute error 0.165525783974 relative to 0.457305783974. It misses both 0.002 and 0.02 tiers. This is a wrong-field/contour error, rather than the larger arithmetic error on a valid problem integral required for the three-point alternative. |
| C — radial integration and thickness convention | 5 / 5 | Correctly integrates W proportional to r^-1 with the area element r dr dtheta, obtaining U(R)=alpha epsilon_0 sigma_0 K^(n+1) Wbar R and C=Wbar/I_n. The per-unit-thickness convention and radial factor are correct independently of the erroneous angular field and J. |
| Verification — independent method or limiting check | 10 / 10 | Unmodified validate3.py checks the n=1 angular ODE against exact biharmonic cosine solutions and reports a maximum absolute difference of 4.44e-16, reproduced on rerun. This is an independently valid analytic limiting check. The final contour comparison shares incorrect parity and does not establish the correctness of H or C. |
| Verification — demonstrated residuals/convergence | 5 / 5 | Unmodified validate.py reproduces upper-half equilibrium residuals 1.45e-5 and 7.24e-5, and a sigma:epsilon versus Q^(n+1) residual of 6.94e-18. These numerical equation checks earn credit; the small sigma_rr(pi) residual does not check the required normal face traction. |

**Total: 55 / 100.** Category subtotals: stress exponent 15/15; mechanics 18/25; H 2/25; C 5/20; verification 15/15.

The reference was rerun from a scratch copy and reproduced its recorded values, collocation/shooting agreement, Cartesian contour checks, and n=1 normalization check. No reference error was found, so an alternative corrected-reference grade is not warranted. The submitted ANSWER.py also reproduces its own reported values; the discrepancies originate in its mechanics and contour evaluation. Rerun outputs are retained in scratch, with detailed diagnostics in [audit_results.json](scratch/audit_results.json).

**(a) Scientific errors**

- At a crack face theta=pi, the face normal is in the theta direction. Free traction requires sigma_theta_theta=sigma_rtheta=0, equivalently phi(pi)=phi'(pi)=0. Enforcing sigma_rr=0 instead leaves the final normalized sigma_theta_theta(pi)=0.718389735, despite sigma_rr(pi)=7.85e-10 and sigma_rtheta(pi)=-1.88e-10.
- Smooth Mode I symmetry requires phi'''(0)=0 as well as phi'(0)=0. The final initial condition phi'''(0)=0.093256796 violates this; the reconstructed normalized U_theta(0)=0.684653117 also violates opening symmetry.
- For u proportional to r^p, the stipulated tensor-shear kinematics give U_r'+(p-1)U_theta=3q sigma_rtheta. The final solver instead uses U_r'-(p+1)U_theta=3q sigma_rtheta, despite the correct later reconciliation paragraph. Evaluating the correct physical shear equation on its normalized field gives a maximum residual of 0.702766104.
- ANSWER.py omits sigma_rtheta times the theta component of u_,x by calling the product odd. Both factors are odd in a Mode I field, so their product is even. robust_j.py instead extends sigma_rtheta evenly, causing an artificial cancellation. On the submitted upper-half field, the omitted full-circle shear-work contribution is 0.473357120; including it changes the contour value from approximately 7.768663168 to 7.295306048. Neither value repairs the field's compatibility and boundary errors.
- The verification's energy identity misses n/(n+1): sigma:epsilon=Q^(n+1), whereas W=n Q^(n+1)/(n+1). The cited 7e-18 residual verifies the former in validate.py. The final energy integration itself includes the correct prefactor, so its constitutive/radial method credit is retained.

**(b) Unsupported claims**

- Agreement of the angular and Cartesian values near I_n=7.7687 does not verify a physical symmetric opening field: both calculations share the incorrect angular solution, and the Cartesian calculation uses even shear-stress parity. Radius independence of a homogeneous ansatz also cannot establish correct boundary conditions or compatibility.
- The claim of H and C agreeing at several radii to nine digits is not demonstrated for a valid final field. The actual Cartesian rerun gives J between 7.768465 and 7.768573 at the four stated radii, while the nine-digit H/C statement has no corresponding successful physical-field output.
- The accompanying indep_bvp.py is not a successful independent solution: its unmodified rerun reports solve_bvp success=False and H=NaN. clean_final.py returns different coefficients (H=2.40743143, C=0.67923413), then fails with an array-shape error. Neither supplies validation of the reported final numbers.

**(c) Presentation problems**

- The main method uses the wrong shear equation, while the reconciliation section identifies the correct equation without updating the final solver or numbers.
- The answer presents conflicting initial curvatures (-0.7050162 versus -0.9142244740), normalized B(0) values (2.827075 versus 1.51569381), and I_n values (7.7687 versus 0.2546) without resolving which angular solution determines the answer.
- Numerous alternative scripts and incompatible normalizations obscure the final method. The definitive ANSWER.py reproduces the requested-number table, but several accompanying checks validate a different expression or fail. These presentation observations do not add deductions outside the rubric.

**Requested-number tolerance check:** p is within 1e-5; H is outside 0.005 (error 0.257350274787); C is outside 0.002 (error 0.165525783974).
