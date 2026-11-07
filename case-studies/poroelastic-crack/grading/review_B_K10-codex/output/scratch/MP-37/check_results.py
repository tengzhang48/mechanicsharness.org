"""Check the delivered schema, printed reproduction, and independent identities."""

import json
from pathlib import Path
import subprocess

import numpy as np

from factor_check import factor


expected = json.loads(Path("results.json").read_text(encoding="ascii"))
completed = subprocess.run(
    expected["reproduce_command"], shell=True, check=True,
    capture_output=True, text=True,
)
printed = json.loads(completed.stdout)
written = json.loads(Path("results.json").read_text(encoding="ascii"))
required = [
    "K1_0", "K1_inf", "V1_0", "V1_inf", "V1_short_coeff", "p_faces_0",
    "K2_0", "K2_inf", "V2_inf",
]
assert expected == printed == written
assert all(isinstance(written[key], (float, int)) for key in required)
assert all(isinstance(value, str) for value in written["formulas"].values())
print("Reproduction: stored JSON == command stdout JSON == regenerated JSON: PASS")
print("Schema: all 9 Report quantities are numbers; all formulas are strings: PASS")

d, eta = .2, 1 / 3
R = written["K1_0"]
energy_A = (
    16 * d / ((1 + eta) * np.sqrt(np.pi))
    + 4 * np.sqrt(np.pi) * d * np.sqrt(2 / (eta * (1 + eta)))
    * np.sqrt(1 - (1 + eta) * R**2)
)
print(f"Volume from independently derived energy identity: {energy_A:.12f}")
print(f"Reported volume coefficient: {written['V1_short_coeff']:.12f}")
print(f"Difference: {energy_A - written['V1_short_coeff']:.3e}")
assert abs(energy_A - written["V1_short_coeff"]) < 1e-7

weak = factor(1e-4, 262144, .01)
weak_expected = 1 - 3e-4 / 4
print(f"Weak coupling eta=0.0001: K={weak['K1_0']:.12f}; first-order K={weak_expected:.12f}")
print(f"Weak-coupling scaled deficit (1-K)/eta: {(1-weak['K1_0'])/1e-4:.9f}; expected limit 0.75")
assert abs((1 - weak["K1_0"]) / 1e-4 - .75) < .001
print("Energy relation and weak-coupling checks: PASS")
