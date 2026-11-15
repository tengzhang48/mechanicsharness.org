# Fracture Mechanics 101 — Agent Edition

**Version:** 0.1.0 (candidate; not yet used in a scored run)  
**Audience:** scientific AI agents with mechanics and numerical tools  
**Purpose:** transferable physical knowledge, not an answer key  
**Scope:** solid mechanics of cracks, including mixed modes, planar curved fronts, and heterogeneous resistance.

This file is a compact reference. It is **not** a prompt to give a particular answer, a verified implementation, or a substitute for the task's boundary conditions. It contains no case-specific reference values, hidden keys, grader outputs, or special crack-front perturbation coefficients.

## A. Definitions and constitutive scope

### A1. Geometry, units, and local basis

- Specify the crack plane, direction of mean advance, and tangent to the crack front. For a non-straight front, define the **local** tangent and in-plane advance direction at each front location.
- Mode I opens the faces; Mode II shears the faces **in the crack plane, perpendicular to the local front**; Mode III shears them **parallel to the local front**. In three dimensions, the labels refer to this local frame, not an arbitrary global coordinate system.
- The crack opening or face slip is the **jump in displacement** across two faces, not the displacement of only one face. State traction, stress, shear and displacement-jump signs.
- Stress intensity K has units Pa√m. Energy release rate G and resistance Γ have units J/m² = N/m. Treat a nominally 2D calculation as per unit out-of-plane thickness when appropriate.
- Stress concentrations may have different asymptotic forms at interfaces, corners, finite process zones, or material nonlinearities. Do not assume an r^(−1/2) singularity without checking applicability.

### A2. Plane stress, plane strain, and bulk material

For an isotropic linearly elastic material of Young's modulus E and Poisson's ratio ν:

- μ = E / [2(1 + ν)] is the shear modulus.
- The effective modulus is E′ = E in plane stress, and E′ = E / (1 − ν²) in plane strain.
- A local segment of an interior three-dimensional crack front in a homogeneous isotropic solid uses the conventional plane-strain singular fields, away from free-surface intersections.
- Distinguish far-field effective properties from material properties near a boundary layer or interface.
- The formulas in this unit are **conditional**; anisotropic, bimaterial, viscoelastic and finite-strain problems require their appropriate constitutive formulations.

## B. Linear elastic crack fields and energy

### B1. Stress intensities

- Under linear elastic fracture mechanics (LEFM), near a smooth crack tip/front and outside a small process zone, stresses scale as K / √(2πr). The angular field and K depend on the loading mode.
- For certain simple geometries, K = Y σ √(πa). The dimensionless geometry factor Y depends on crack shape, boundary conditions, location and load; do not assume Y = 1 outside its idealized domain.
- Crack opening typically scales as (K / E′)√r; the coefficient depends on the precise geometry and definition of the displacement jump.
- A uniform far-field stress, applied crack-face pressure, a face-bridging traction and an internal pore pressure are different loadings. Use their correct boundary conditions.

### B2. Energy release rate

For a local segment of a planar crack in a homogeneous isotropic elastic three-dimensional body, under the standard singular-field conventions:

    G = (K_I² + K_II²) / E′ + K_III² / (2μ)

Here E′ = E/(1 − ν²) for the local interior plane-strain field. Pure plane-stress, two-dimensional modeling uses its corresponding effective modulus and admissible modes. This simple diagonal quadratic form is not universal across material interfaces or anisotropic media.

- G is elastic energy available per increment of crack area; Γ is the resistance required to create that area.
- A local condition G = Γ can characterize quasistatic advance under suitable assumptions, but neither the equality nor a computed equilibrium alone guarantees stability.
- The appropriate energy expression depends on load control versus displacement control, crack-tip regularization, and any bulk dissipation.
- For a conservative elastic body, configurational force is obtained from the variation of the appropriate total potential energy with crack geometry. This underpins reciprocity and some independent consistency tests.

### B3. Superposition and crack-face loads

- In linear elasticity, superposition is valid for stresses, displacements and stress intensity factors at fixed geometry and material properties.
- It does **not** generally hold for nonlinear material laws, unilateral contact, or evolving cohesive damage without additional reasoning.
- A closure or bridging traction may oppose applied opening. Distinguish its sign on the faces from its positive work of separation.
- Weight functions can map distributed crack-face loads to tip fields when the governing linear boundary-value problem supports them. A remembered weight function must not be transplanted to a different geometry without justification.

## C. Fracture mechanisms beyond LEFM

### C1. Cohesive and bridging zones

- A traction–separation law relates resistance to the full face-to-face gap or slip; its area is the work of separation per area.
- The **intrinsic tip resistance** and work consumed in a wake or bridging zone are not automatically interchangeable.
- Small-scale bridging may admit a macroscopic effective toughness when the zone is short relative to relevant dimensions. For large-scale bridging, solve the elastic body and the bridging law together.
- A material length built from modulus, strength and separation is a scaling estimate, not necessarily the actual process-zone length.
- Zero singular stress intensity at an artificial outer yield-zone tip does not imply zero opening elsewhere in the crack.

### C2. Plasticity, J, and HRR-type fields

- Small-scale yielding permits a far-field elastic K description even if the immediate crack tip yields.
- For appropriate nonlinear elastic or deformation-theory settings, J is the energy flux/configurational integral; its path independence requires its standard hypotheses, including absence of unaccounted body forces, inelastic history or other sources inside contours.
- A power-law hardening material can have near-tip HRR-type stresses proportional to r^(−1/(n+1)) when strain scales as stress to the power n; the corresponding strain exponent is −n/(n+1).
- The dimensionless angular functions, material normalization, amplitude and any eigenvalues must come from the specific boundary-value problem. The scaling exponent alone is not a complete solution.
- Check whether a model assumes plane stress, plane strain, monotonic deformation, small strains or a particular hardening convention.

### C3. Poroelastic cracks and diffusion

- A fluid-saturated porous solid couples solid deformation, pore pressure and fluid flux. Applied face pressure is not identical to pore pressure within the skeleton.
- Under fast loading, large regions may behave undrained; under slow loading, drained constitutive response may apply. These are limiting conditions, not two independently assigned materials.
- A diffusion length grows as √(ct), with diffusivity c. The comparison between this length, the crack size and the distance from a drained boundary matters; a single time label need not describe every spatial scale.
- A singular-tip observable may require a separate asymptotic analysis from crack-scale displacement or volume. The correct near-tip and time limits depend on the actual boundary conditions and sequence of limits.
- Check flux and pressure boundary conditions, mass balance and dimensional consistency.

### C4. Viscoelastic fracture

- Viscoelastic stress depends on deformation history (or strain depends on stress history); instantaneous and relaxed moduli need not be equal.
- Creep compliance and stress-relaxation modulus are distinct response functions. Their characteristic times need not coincide.
- A steadily growing crack can dissipate energy in the bulk while a cohesive/process zone consumes intrinsic fracture work. A far-field energy release rate may include both.
- Check the assumed steady traveling frame, admissible propagation speeds, material thermodynamics and finite-geometry constraints. A high-speed asymptote need not be attained in a finite strip.

### C5. Films and diffusion wedges

- Prestress is part of the **total traction** that must vanish on a traction-free crack face; solving for only the stress change can produce an incorrect face condition.
- A thin film bonded to a substrate may permit a reduced membrane/spring-bed model, but its length scale and energy release depend on that reduction's assumptions.
- Diffusion-assisted opening of a grain boundary couples chemical potential, transport and elasticity. Long-time crack-like behavior and relaxation time require a model-specific eigenproblem or asymptotic limit.
- A time scale formed by dimensional analysis need not equal the slowest actual e-folding time. Never infer numerical eigenvalues from scaling alone.

## D. Mixed-mode three-dimensional crack fronts

### D1. Distinguish two sources of response to front shape

For a planar front that varies along its length, even a small perturbation can change local stress intensity factors through:

1. **Nonlocal elastic redistribution:** a change in front position at one location changes the field at other locations.
2. **Local frame/orientation effects:** the local tangent and advance direction change; the local Mode II/III decomposition must reflect the actual three-dimensional elastic solution.

These contributions are conceptually distinct. A simple rotation of the original K vector is not, by itself, a full first-order solution. Nor is copying the Mode I perturbation kernel into both shear modes generally valid.

### D2. Why Modes II and III are not interchangeable

- Shear slip aligned with its spatial variation (edge-like) and transverse shear slip (screw-like) probe different elastic compliances in a three-dimensional isotropic body.
- Curvature or waviness can mix the shear modes and expose this distinction, even if the reference straight front carries a simple prescribed mixture.
- The response can depend on ν, geometry, load, and the actual local frame. There is **no universal scalar line-tension coefficient** that replaces the three-dimensional elasticity problem.
- Establish the domain of perturbation theory: nearly planar crack, small displacement and slope relative to relevant scales, homogeneous medium if assumed, and a well-defined straight reference front.
- For larger changes, front kinks, cracks leaving the plane or nonlinear process zones, reassess the model rather than extrapolating first-order formulas.

This reference intentionally gives **no** special mixed-mode crack-front coefficient formulas, tabulated kernels, or solutions to a spatially prescribed front.

### D3. Trapping by heterogeneous resistance

- Let G(s; a, loading) be the elastic driving force along a front with shape a; let Γ(s, state) be local resistance.
- Where the front advances quasistatically under a local rate-independent law, G = Γ may be required. Pinned portions can satisfy a no-growth inequality G ≤ Γ, depending on the stated irreversibility law.
- Because front segments interact elastically, prescribing Γ(s) does not determine the shape without solving how a changes G.
- Distinguish existence of an equilibrium front from its stability and reachability under a loading history.
- The phrase "line tension" is a **qualitative analogy**, not automatically a quantitatively accurate local crack-front law. A nonlocal first-order response is different from a universal short-range tension.

## E. Applicability boundaries

Do not silently apply these textbook relations to:

- Oscillatory elastic interface singularities, anisotropic media or graded materials without the corresponding tip theory.
- Finite deformation, nonlinear elasticity or plastic loading/unloading without the correct constitutive framework.
- Fast dynamic fracture, inertia-dominated fields, or rate-dependent fracture laws using quasistatic G as though no extra terms existed.
- Fronts undergoing strong deflections, self-contact, kinking or changing fracture planes using small-slope planar perturbation formulas.
- A crack with a process zone comparable to the problem dimensions while assuming small-scale fracture.

For any particular scientific task, first identify the actual governing equations, boundary conditions, constitutive model, geometric scales, and target observables. Any special coefficient, eigenvalue or profile not stated here must be derived, computed or otherwise justified within the allowed tools. Do not assert an unperformed calculation as verification.

---

**Content boundary:** This physics reference contains no benchmark tasks or grading keys and is intentionally less specific than the human study pages. It is designed for controlled evaluation in an isolated agent workspace. A separate file, `verification-workflow.md`, describes an optional explicit verification intervention and should not be silently bundled with this document when comparing physics-only against physics-plus-verification conditions.
