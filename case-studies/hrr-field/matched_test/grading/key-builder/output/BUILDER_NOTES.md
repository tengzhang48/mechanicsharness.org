# Builder notes

This set contains three fresh probes. Only each folder's `TASK_PROBLEM.md` is intended for the solver. `REFERENCE.md`, `solve_reference.py`, `reference_output.json`, and `RUBRIC.md` are grading/build artifacts. All calculations use mm, MPa, and N, with explicit tensor/engineering-shear conventions. The reference code needs only numpy and scipy (A1 needs numpy only).

## V1 — out-of-plane shear, slip, and strain contours

V1 changes the loading mode and governing field from an in-plane tensile system to a scalar antiplane displacement problem. The specified homogeneous shear law is the local power-law mechanism used in the source, adapted to engineering shear. The solver must formulate equilibrium and compatibility, infer the admissible radial behavior, normalize the field by energy flux, and calculate face slip, signed off-axis stresses, and strain-contour geometry.

This is not a repeat of the earlier plane-strain n = 6.5 stress exponent/amplitude/energy-ratio probe. The requested outputs are displacement, two shear components, and strain contours; the angular problem is second order rather than the in-plane fourth-order stress-function problem. It is also not a notch-angle variation. The source itself discusses antiplane shear as a neighboring problem when assessing the scope of local similarity solutions; it does not supply this problem's angular solution or numerical data.

The pure shear law is deliberately specified without an elastic term, so the task does not hide an elastic-to-plastic matching calculation. The given leading-field dominance and energy flux make the local problem complete. Approximate expected effort: 10–15 minutes to formulate, 10–15 to solve the angular problem/normalize, and 10 to calculate and verify. Both a numerical IVP/BVP route and an exact phase parametrization are valid. A capable solver need not know the exact phase solution.

Useful grader distinctions: engineering shear is ∇w; W is n/(n+1) times τ·γ; face traction vanishes while radial shear need not. The initiation question requires a failure law but does not justify refusing the stationary calculation. All seven numerical references are independently checked by the angular ODE and a phase parametrization/exact flux, with a further n = 1 linear special-solution check.

## V2 — inverse calibration with a different material asymptote and constraint

V2 uses the source's alternative of a positive finite high-stress tangent stiffness, rather than a power-law compliance, and uses a plane-stress membrane. It turns the forward local-field analysis into an inverse problem: infer the unknown tensile slope from the transverse strain coefficient, then predict other strain/displacement coefficients and a screened radial interval.

The central transfer is to preserve the elastic bulk response while the deviatoric compliance changes. The effective leading isotropic pair is (T, ν_t), with ν_t tending toward 1/2 as T/E becomes small. Holding the original ν fixed while replacing E by T gives a different answer. The bounded strain contribution from the yield intercept does not change the requested leading coefficients. The supplied full tensor law removes any ambiguity about how a uniaxial bilinear curve is generalized to multiple stress components.

This is distinct from all three earlier probes: no hardening exponent is to be fitted or reported, there is no notch, and the outputs emphasize inverse material identification, transverse strain, compatible displacement, and an explicit nonempty radius screen. Expected effort: roughly 15 minutes for the tensor reduction and local field, 5 for inversion, and 15–20 for outputs and a check.

The coefficient measurement is explicitly a formal small-strain computational limit. The radial screens are explicitly leading-field tests. Their numerical tolerances are **not** physical error bounds for the complete finite-r solution, especially at the high-stress screen where σ_y/σ_e = 0.1. This distinction avoids asking for an unprovided global nonlinear boundary-value solution. There is no supplied failure criterion. The references are checked by a tensor-based inverse/root calculation, direct flux integration, displacement/strain compatibility, uniaxial calibration, and the T = E limit.

## A1 — mechanism, volume change, and three-dimensional constraint

The three decisions are: applicability after antiplane adaptation; failure of a traceless-plastic-strain law to describe cavity-induced volume change; and uncertainty in a single 2D constraint near a finite-thickness crack front. The set tests both avoiding unnecessary refusal and recognizing missing physics or missing validation. It does not repeat the earlier cyclic-loading, large-scale-yielding, or contour-agreement applicability set.

Case 1 has favorable length scales and quantitatively modest principal strains with dominant plastic shear. Its model applies even though the loading is not tensile. Case 2 contains a direct constitutive contradiction: elastic bulk strain is 0.006 while observed volume strain is 0.018, requiring 0.012 plastic dilation. Cavity/porosity evolution and pressure-sensitive flow are relevant; an extra fracture criterion alone cannot repair the requested volume response. Case 3 preserves the material mechanism but has r/(distance to surface) = 0.30–1.00 and yielded extent/thickness = 1.25. It calls for constraint validation or a 3D field rather than automatic rejection of the material law.

A1's reference code checks tensor versus scalar shear and tensor versus bulk-modulus volume response. It audits geometry in two unit systems and checks the distinct plane-stress/plane-strain conditions in a known elastic limit. No actual 3D mid-plane stress or 15% approximation error has been calculated; the problem intentionally lacks enough information to certify that accuracy. The reference states that limitation explicitly. Equivalent conditional wording must receive full decision credit. False applicability and unnecessary refusal carry equal decision-point penalties.

## Relationship to the supplied source and build audit

The local power-law mechanism, finite-work requirement, and energy-flux normalization are discussed in the source's singular-field and amplitude sections (images `pages/page-04.png` through `pages/page-11.png`). The alternative piecewise-linear material and its plane-stress reduction are in `pages/page-12.png` and `pages/page-13.png`. The source's assumptions about proportionality, elastic bulk response, small kinematics, and constraint motivate A1 (especially `pages/page-02.png`, `pages/page-16.png`, and `pages/page-17.png`). Equations in the statements are fully dimensional and independently defined; they do not point solvers to source authors or equation numbers.

References were computed locally with newly written code, without web search/fetch or outside resources. Every numeric variant value has a second-route check and a known special/material limit. Case 3's intentionally unavailable 3D accuracy is expressly marked as uncomputed. Each rubric totals 100 points and gives explicit partial credit, including consequential-error credit. The time-stamped steps are in `PROCESS.md`; all build files and scratch work remain inside this folder.
