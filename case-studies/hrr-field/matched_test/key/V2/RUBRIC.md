# Rubric — 100 points

Use the final leading-field quantities as the scoring target. A full finite-r boundary-value solution is not required. Any correct analytical, numerical, or variational method earns full credit; reference-code reuse is not required.

| Category | Points |
|---|---:|
| Model choice and justification | 30 |
| Numerical correctness | 40 |
| Completion and conventions | 10 |
| Quantitative verification | 15 |
| Approximation scope and missing physics | 5 |
| **Total** | **100** |

## Model choice and justification — 30

- **10:** Reduces the high-stress leading tensor compliance correctly, preserving the elastic bulk response and allowing nonzero ε_zz. Full credit for the equivalent isotropic pair (T, ν_t), with ν_t = 1/2 − (1/2−ν)T/E. Give 5 for correctly retaining plastic incompressibility/elastic volume response but making one algebraic error; 2 for recognizing that T is relevant while incorrectly holding ν fixed; 0 for imposing ε_zz = 0 contrary to the given constraint.
- **10:** Chooses and justifies the appropriate traction-free tensile leading field of the reduced law. Give 5 for the correct stress/displacement radial behavior with an incomplete explanation or angular error; 0 for an unrelated power-law exponent or field. A known linear crack solution is a valid method after the constitutive reduction is established.
- **5:** Correctly normalizes the local field using J = K_loc²/T, or an equivalent direct flux calculation. Give 3 for the right flux argument with one factor error; 0 for using only E without adapting the local field.
- **5:** Uses the transverse measurement to solve for an admissible unique T and relates ε_zz to the in-plane stresses. Give 3 for the correct inverse equation without completing the inversion, and 0 for treating C_zz as a direct axial compliance or ignoring the measurement.

## Numerical correctness — 40

| Output | Points | Full-credit reference | Tolerance |
|---|---:|---:|---:|
| T | 12 | 4816.472187 MPa | 1% |
| C_rr at 60° | 5 | +0.003942048238 mm^1/2 | 2% |
| C_θθ at 60° | 5 | +0.0006331781351 mm^1/2 | 2% |
| C_γrθ at 60° | 5 | +0.005731131133 mm^1/2 | 2% |
| C_y at 90° | 5 | +0.009139018941 mm^1/2 | 2% |
| Lower radius endpoint | 4 | 0.01000000000 mm | 2% |
| Upper radius endpoint | 4 | 0.06132522855 mm | 2% |

Apply relative tolerances. An outside-tolerance value earns half its points if the relevant coefficient or screening condition is correctly computed using a single identifiable upstream error; otherwise give 0. Correct formulas without numerical evaluation earn half the corresponding points. The lower screen must examine all three principal strains, including ε_zz; using only in-plane principal strains is an additional error. Wrong signs and tensor/engineering shear confusion cannot earn full numerical credit. The interval must be explicitly recognized as nonempty to earn full endpoint credit.

## Completion and conventions — 10

- **5:** All seven requested numerical outputs are reported and identifiable; deduct 1 per absent output, down to 0.
- **5:** Correct units (2), tensor/engineering shear convention (2), and consistent polar/displacement signs (1).

## Verification — 15

Award 5 for a genuine independent condition, 5 for its quantitative evaluation, and 5 for interpreting the agreement. One well-executed check can receive all 15: for example, J from W and the displacement field; sym(∇u) versus the original leading tensor law; a separate tensor calibration; or the T = E known limit. Substitution back into the measured C_zz alone earns at most 3 because it merely reuses the calibration equation. A dimensions-only check earns at most 5.

## Scope and missing physics — 5

- **3:** Recognizes that the three screens give a working interval for the stated leading approximation, without proving full-solution error or physical validity at every smaller radius. Give 1 for identifying only a small-strain or continuum limitation; 0 for equating the 2% numerical tolerance with physical finite-r accuracy.
- **2:** States that initiation requires an additional failure/fracture/damage criterion. Give 1 for saying that initiation is undetermined without explaining why; 0 for inferring a failure load from the strain screen alone.

An unjustified refusal because the formal field diverges loses credit for the work not done: the problem explicitly asks for coefficients and a screened interval. This does not prevent a solver from explaining genuine limits after completing the calculation.
