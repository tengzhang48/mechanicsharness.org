"""Read-only numerical audit of the copied submission and grading key."""
import contextlib
import io
import json
import runpy

import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    sub = runpy.run_path("solve.py")

with open("reference_rerun.json") as stream:
    key = json.load(stream)["route_B_phase_formula"]

reported = [
    ("slip_0.02_mm", 0.0013181),
    ("slip_0.08_mm", 0.0017937),
    ("tau_rz_0.02mm_60deg_MPa", 38.59),
    ("tau_thetaz_0.02mm_60deg_MPa", 138.01),
    ("gamma_contour_0deg_mm", 0.04689),
    ("gamma_contour_90deg_mm", 0.02507),
    ("gamma_contour_180deg_mm", 0.01340),
]
print("Reported numbers versus independently rerun reference:")
for name, value in reported:
    err = abs(value - key[name]) / abs(key[name])
    print(f"{name}: answer={value:.12g}, key={key[name]:.12g}, error={100*err:.9f}%")
    assert err <= 0.02

print("\nActual stress implementation versus displacement compatibility:")
r = 0.02
h = 1e-7
for degrees in (60.0, -60.0):
    theta = np.deg2rad(degrees)
    x, y = r * np.cos(theta), r * np.sin(theta)
    grad = np.array([
        (sub["w_of"](x+h, y) - sub["w_of"](x-h, y)) / (2*h),
        (sub["w_of"](x, y+h) - sub["w_of"](x, y-h)) / (2*h),
    ])
    expected = 80 * (np.linalg.norm(grad) / 0.002)**(1/3.5) * grad / np.linalg.norm(grad)
    actual = np.array(sub["tau_cart"](x, y))
    print(f"theta={degrees:+g} deg: w={sub['w_of'](x,y):+.12g} mm")
    print(f"  grad w={grad}; stress implied by grad w={expected} MPa")
    print(f"  tau_cart={actual} MPa; relative vector mismatch={np.linalg.norm(actual-expected)/np.linalg.norm(expected):.9g}")

axis_tau = sub["tau_cart"](r, 0.0)
expected_axis = 80 * (sub["field"](r, 0.0, sub["A_amp"])[3] / 0.002)**(1/3.5)
print(f"theta=0: tau_cart={axis_tau} MPa; expected tau_yz={expected_axis:.12g} MPa")

state60 = sub["field"](r, np.pi/3, sub["A_amp"])
print(f"\nEngineering shear at requested stress point: {state60[3]:.12g} ({100*state60[3]:.8f}%)")
print(f"F(pi)={float(sub['F_s'](np.pi)):.12g}")

upper_j = sub["J_circle"](r, sub["A_amp"])
print(f"Submission's numerical flux={upper_j:.16g} N/mm")
print("The lower-half sign error reverses both stress and assumed gradient, preserving W and tau_r*w_x in this integral.")
