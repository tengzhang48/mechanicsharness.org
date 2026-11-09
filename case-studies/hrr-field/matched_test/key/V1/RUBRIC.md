# Rubric — 100 points

Score any correct formulation; reference-code reuse, a particular normalization, and a particular solver are not required. Apply tolerances to the requested signed quantities. Award consequential-error credit as specified below rather than repeatedly zeroing correct downstream work.

| Category | Points |
|---|---:|
| Model choice and physical justification | 30 |
| Numerical correctness | 40 |
| Completion and conventions | 10 |
| Quantitative verification | 15 |
| Recognition of missing fracture physics | 5 |
| **Total** | **100** |

## Model choice and justification — 30

- **8:** Uses antiplane displacement compatibility, equilibrium, and the supplied nonlinear shear law. Full credit for an equivalent stress or variational formulation. Give 4 for the correct nonlinear mechanism with an incomplete governing formulation; 0 for substituting an unrelated tensile field or an unqualified linear shear law.
- **8:** Obtains the admissible local radial behavior from a nonzero finite energy flux, the nonlinear eigenproblem, or an equivalent valid argument. Give 4 for a plausible separable nonlinear field with an unsupported radial exponent; 0 for an incompatible exponent. The exponent need not be separately reported if the correct reasoning/field is clear.
- **6:** Implements odd displacement symmetry, the traction-free face condition, and the positive-slip loading convention. Give 2 per item; partial credit for a clearly correct intent with a transcription error.
- **8:** Fixes field magnitude from the supplied J with the correct work density and stress/displacement contribution. Give 4 for the correct flux formulation with one normalization error, 2 for identifying the load quantity but not using it, and 0 for leaving magnitude arbitrary. Confusing work with complementary work does not earn full credit.

## Numerical correctness — 40

| Output | Points | Full-credit reference |
|---|---:|---:|
| Δw at 0.0200 mm | 8 | +0.001318126778 mm |
| Δw at 0.0800 mm | 8 | +0.001793693738 mm |
| τ_rz at 60°, 0.0200 mm | 6 | +38.58659686 MPa |
| τ_θz at the same point | 6 | +138.0106166 MPa |
| Contour radius at 0° | 4 | 0.04689442305 mm |
| Contour radius at 90° | 4 | 0.02506612349 mm |
| Contour radius at 180° | 4 | 0.01339840659 mm |

Each number earns full credit within 2% relative error. For a value outside tolerance, give half its points if the relevant quantity is correctly defined and evaluated consistently from a field with one identifiable upstream error; otherwise give 0. A correct derivation/formula without numerical evaluation earns half the corresponding points. A wrong sign cannot earn full credit. Do not award numerical credit for unexplained guesses outside tolerance.

## Completion and conventions — 10

- **5:** All seven numerical outputs are present and identifiable; deduct 1 per absent output, down to 0.
- **5:** Correct mm/MPa units (2), engineering versus tensor shear (2), and explicit face-limit convention (1). Recognizing the face limit in the calculation suffices; an extra sentence is not mandatory.

## Verification — 15

One sound quantitative check can earn all 15 points. Award 5 for a check that addresses a genuine field, boundary, energy, or limiting-condition requirement; 5 for showing a numerical result or residual; and 5 for interpreting its agreement with the required value/limit. Examples include J on a contour, traction residual plus symmetry, the independently known linear limit, or an independent solution method. A restatement of input J or a claim that software converged earns no more than 3. Checking only dimensions earns at most 5.

## Missing physics — 5

Full credit for saying that an initiation load needs fracture resistance, a damage law, a cohesive/critical-separation rule, or an equivalent failure criterion. Give 2 for saying initiation is undetermined without identifying the missing ingredient; 0 for claiming that the stationary constitutive field alone fixes initiation. A refusal to perform the supplied stationary calculation receives no credit in the uncompleted categories.
