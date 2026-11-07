"""Independent grading calculations; does not modify the submission."""

import cmath
import json
import math
from pathlib import Path

root = Path(__file__).resolve().parents[1]
original = json.loads((root / "submission/results.json").read_text())
rerun = json.loads((root / "scratch/reproduction/results.json").read_text())
nu, nu_u, alpha = 0.2, 0.4, 0.89
eta = alpha * (1 - 2 * nu) / (2 * (1 - nu))
reference = {
    "K1_0": 0.807593,
    "K1_inf": 1.0,
    "V1_0": math.pi * (1 - nu_u),
    "V1_inf": math.pi * (1 - nu),
    "V1_short_coeff": 2.440220,
    "p_faces_0": 2 * (nu_u - nu) / (alpha * (1 - 2 * nu)),
    "K2_0": 0.0,
    "K2_inf": -eta,
    "V2_inf": -alpha * math.pi * (1 - 2 * nu) / 2,
}

print("Original JSON equals rerun JSON:", original == rerun)
for name, expected in reference.items():
    claimed = original[name]
    tolerance = 0.01 if expected == 0 else 0.01 * abs(expected)
    print(
        f"{name}: submission={claimed:.12g}, key={expected:.12g}, "
        f"pass={abs(claimed - expected) <= tolerance}"
    )

face = 16 * (1 - nu_u) * (nu_u - nu) / ((1 - nu) * math.sqrt(math.pi))
print(f"Face-only reference A: {face:.12g}")
print(f"Submitted A / face-only reference A: {original['V1_short_coeff'] / face:.12g}")
print(f"Combined-loading K_inf: {1 - eta:.12g}")
print(f"Combined-loading V_inf: {reference['V1_inf'] + reference['V2_inf']:.12g}")
print("Alternative effective-stress loading is alpha*loading_1 + loading_2:")
print(f"K_0: {alpha * reference['K1_0']:.12g}")
print(f"K_inf: {alpha - eta:.12g}")
print(f"V_inf: {alpha * reference['V1_inf'] + reference['V2_inf']:.12g}")

print("Harmonic pressure with p=1 on slit and p=0 on a confocal outer ellipse:")
xi = cmath.acosh(2 + 0j).real
for xi_R in [2, 4, 8, 16, 32]:
    print(f"xi_R={xi_R}, p_at_x_2L={1 - xi / xi_R:.9f}")

print("Key COD-to-K dimensional typo; L=1 and unit opening traction:")
for G in [1.0, 4.0]:
    cod_coefficient = 2 * (1 - nu) * math.sqrt(2) / G
    printed = math.sqrt(math.pi * G) * cod_coefficient / (math.sqrt(8) * (1 - nu))
    corrected = G * math.sqrt(math.pi) * cod_coefficient / (math.sqrt(8) * (1 - nu))
    print(
        f"G={G:g}, printed_key_K={printed:.12g}, "
        f"corrected_K={corrected:.12g}, expected={math.sqrt(math.pi):.12g}"
    )
