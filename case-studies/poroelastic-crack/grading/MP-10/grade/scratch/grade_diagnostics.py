"""Read-only numerical diagnostics of the delivered formulas and boundary conditions."""
import cmath
import json
import math
from pathlib import Path

from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
delivered = json.loads((ROOT / "submission/results.json").read_text())
rerun = json.loads((ROOT / "scratch/submission/results.json").read_text())
nu, nu_u, alpha, L, G = 0.2, 0.4, 0.89, 1.0, 1.0
reference = {
    "K1_0": 0.807593,
    "K1_inf": 1.0,
    "V1_0": math.pi * (1 - nu_u),
    "V1_inf": math.pi * (1 - nu),
    "V1_short_coeff": 2.440220,
    "p_faces_0": 2 * (nu_u - nu) / (alpha * (1 - 2 * nu)),
    "K2_0": 0.0,
    "K2_inf": -alpha * (1 - 2 * nu) / (2 * (1 - nu)),
    "V2_inf": -alpha * math.pi * (1 - 2 * nu) / 2,
}

print("Reproduced JSON identical:", delivered == rerun)
for key, target in reference.items():
    actual = delivered[key]
    passes = abs(actual - target) <= (0.01 if target == 0 else 0.01 * abs(target))
    print(f"{key}: delivered={actual:.12g}, key={target:.12g}, within band={passes}")

def fp(z):
    return (2 / math.pi) * cmath.asin(L / z)

print("\nLoading-2 pressure on both halves and the crack center (upper/lower):")
for x in [-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9]:
    print(f"x={x:+.1f}: p_up={fp(complex(x, 1e-9)).real:+.9f}, "
          f"p_down={fp(complex(x, -1e-9)).real:+.9f}, prescribed=+1")

def contour_integrand(theta):
    z = 2 * L * cmath.exp(1j * theta)
    return fp(z) * 1j * z

loop_integral = complex(
    quad(lambda theta: contour_integrand(theta).real, 0, 2 * math.pi)[0],
    quad(lambda theta: contour_integrand(theta).imag, 0, 2 * math.pi)[0],
)
swelling_factor = alpha * (1 - 2 * nu) / (2 * G)
print("\nIntegral of f'(z) dz around the slit:", loop_integral)
print("Expected residue value 4iL:", 4j * L)
print("Displacement change after one circuit:", swelling_factor * loop_integral)
print("Thus no single-valued stress-free displacement u=swelling_factor*f exists.")

face_imag_integral = -4 / math.pi * quad(lambda x: math.acosh(L / x), 0, L)[0]
print("Upper-face u_y(L)-u_y(-L):", swelling_factor * face_imag_integral)
print("This cannot be zero at both tips by adding a constant to f.")

def zw(z):
    return -(1 - z / cmath.sqrt(z * z - L * L))

print("\nWestergaard code's principal-square-root branch at far left/right:")
for x in [-100.0, 100.0]:
    print(f"x={x:+g}: Re Z={zw(complex(x, 1e-9)).real:+.9f}, required far-field=0")

faces = 16 * (1 - nu_u) * (nu_u - nu) / (math.sqrt(math.pi) * (1 - nu))
print("\nKey face-only coefficient:", faces)
print("Key tip coefficient:", reference["V1_short_coeff"] - faces)
print("Delivered coefficient / key face-only coefficient:", delivered["V1_short_coeff"] / faces)
print("Claimed (ct/L^2)^(3/4) divided by the leading (ct/L^2)^(1/2) scales as (ct/L^2)^(1/4), not (ct/L^2)^(1/2).")
