"""Generate results.json from the analytical energies and converged FE phases."""
import cmath
import json
import math

# Non-dimensional units: E2=h=1, H=2.  See method.md for the derivations.
EPS = math.log(2.0 / 3.0) / (2.0 * math.pi)
GA = 287.0 / 162.0
GB = 34.0 / 27.0
KABS_A = math.sqrt(1435.0 / 384.0)
KABS_B = math.sqrt(85.0 / 32.0)

# Central values of the converged Williams-series fits.  The opening and
# interface-traction fits bracketing these choices are listed in method.md.
PSI_A = 52.15
PSI_B = 4.24

k_sum = (KABS_A * cmath.exp(1j * math.radians(PSI_A)) +
         KABS_B * cmath.exp(1j * math.radians(PSI_B)))
result = {
    "eps": EPS,
    "caseA": {
        "G_E2_h_over_P2": GA,
        "Kabs_sqrt_h_over_P": KABS_A,
        "psi_deg": PSI_A,
    },
    "caseB": {
        "G_E2_h3_over_M2": GB,
        "Kabs_h32_over_M": KABS_B,
        "psi_deg": PSI_B,
    },
    "psi_both_deg": math.degrees(cmath.phase(k_sum)),
}

with open("results.json", "w") as f:
    json.dump(result, f, indent=2)
    f.write("\n")
print(json.dumps(result, indent=2))
