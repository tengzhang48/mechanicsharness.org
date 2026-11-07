"""Local grading evidence; uses the supplied corrected key, no external sources."""
from pathlib import Path
import json
import math

root = Path(__file__).resolve().parents[1]
reported = json.loads((root / "submission/results.json").read_text())
nu, nuu, alpha = 0.2, 0.4, 0.89
eta = alpha * (1 - 2 * nu) / (2 * (1 - nu))
expected = {
    "K1_0": 0.807593,
    "K1_inf": 1.0,
    "V1_0": math.pi * (1 - nuu),
    "V1_inf": math.pi * (1 - nu),
    "V1_short_coeff": 2.440220,
    "p_faces_0": 2 * (nuu - nu) / (alpha * (1 - 2 * nu)),
    "K2_0": 0.0,
    "K2_inf": -eta,
    "V2_inf": -alpha * math.pi * (1 - 2 * nu) / 2,
}
print("Against the corrected key (1% nonzero; absolute 0.01 for K2_0):")
for name, want in expected.items():
    got = reported[name]
    error = abs(got - want)
    tolerance = 0.01 * abs(want) if want else 0.01
    relative = error / abs(want) if want else error
    print(f"{name}: got={got:.12g}, key={want:.12g}, error={relative:.6g}, passes={error <= tolerance}")

face_paper = 16 * (1 - nuu) * (nuu - nu) / ((1 - nu) * math.sqrt(math.pi))
local_with_correct_pressure = 8 * (nuu - nu) / ((1 - nu) * math.sqrt(math.pi))
print(f"\nCorrect pressure, still using the submission's local opening rule: A={local_with_correct_pressure:.12g}")
print(f"Key's face-only reciprocal contribution: A_face={face_paper:.12g}")
print(f"Key's additional tip contribution: A_tip={2.440220 - face_paper:.12g}")

# An elastic counterexample to 'pointwise uniform pressure implies zero stress'.
# Prescribe p=P in a disk r<R and p=0 outside. In plane strain let u=a*r
# inside, u=b/r outside. Displacement and radial-traction continuity give
# b=a*R^2 and 2G*a/(1-2nu)-alpha*P = -2G*a.
# Hence a=alpha*P*(1-2nu)/(4G*(1-nu)), sigma_inside=-eta*P,
# independent of R. Every finite-R solution has zero stress at infinity;
# p tends to P pointwise as R grows, but the inner stress stays nonzero.
G, P = 1.0, 1.0
a = alpha * P * (1 - 2 * nu) / (4 * G * (1 - nu))
inside_stress = 2 * G * a / (1 - 2 * nu) - alpha * P
print("\nExpanding pressurized disk: an equilibrium limit, not a transient crack simulation")
for radius in (10.0, 100.0, 1000.0):
    b = a * radius**2
    displacement_jump = a * radius - b / radius
    outside_radial_stress = -2 * G * b / radius**2
    traction_jump = inside_stress - outside_radial_stress
    print(f"R={radius:g}: inner stress={inside_stress:.12g}, displacement jump={displacement_jump:.3g}, traction jump={traction_jump:.3g}")
print(f"Crack-face correction in this limiting compressive field: K=-eta={-eta:.12g}, V=-pi*(1-nu)*eta={-math.pi*(1-nu)*eta:.12g}")
