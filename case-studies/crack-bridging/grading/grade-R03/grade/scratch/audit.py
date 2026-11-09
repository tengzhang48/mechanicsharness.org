"""Local grading audit; all computations and artifacts stay in the task folder."""

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent / "submission"))
import explore
import reference

ep, a, t, dc, g0 = explore.Ep, explore.a, explore.t, explore.dc, explore.G0
k0 = math.sqrt(ep * g0)


def stress(c):
    return 2 * t / math.pi * math.acos(a / c) + k0 / math.sqrt(math.pi * c)


def correct_opening(c, sigma):
    return (
        4 * math.sqrt(c * c - a * a)
        * (sigma - 2 * t / math.pi * math.acos(a / c))
        + 8 * t * a / math.pi * math.log(c / a)
    ) / ep


length = brentq(
    lambda ell: explore.delta_a_closed(a + ell, stress(a + ell)) - dc,
    0.0, 1.0, xtol=1e-14,
)
c = a + length
sigma = stress(c)
submitted_sigma, submitted_length = 158.525, 34.324
target_sigma, target_length = 133.730286053811, 4.08560913824553

k = k0 / (t * math.sqrt(math.pi * a))


def slope(u):
    return 2 / (math.pi * math.sqrt(math.expm1(2 * u))) - 0.5 * k * math.exp(-u / 2)


u_peak = brentq(slope, 1.0, 10.0)
u_wrong = math.log(c / a)
u_reference = reference.analytic_solution()["log_eta"]
u_sample = np.linspace(1e-8, u_wrong, 10001)

example_c, example_sigma = a + 1e-3, 100e6
K_exact = explore.K_closed(example_c, example_sigma)
K_quad = explore.K_num(example_c, example_sigma)

# An additional contradiction in the submitted purported full opening field:
# at x=0 it predicts a negative gap for this opening-load example.
submitted_center_gap = 4 / ep * (
    (example_sigma - t) * example_c + t * (0 - example_c + example_c)
)

result = {
    "reported_number_errors_percent": {
        "stress_precise": 100 * (submitted_sigma / target_sigma - 1),
        "length_precise": 100 * (submitted_length / target_length - 1),
        "stress_bold": 100 * (158.5 / target_sigma - 1),
        "length_bold": 100 * (34.3 / target_length - 1),
    },
    "acceptable_intervals": {
        "Sigma_MPa": [target_sigma * 0.997, target_sigma * 1.003],
        "L_mm_per_end": [target_length * 0.997, target_length * 1.003],
    },
    "submission_code_root": {
        "Sigma_MPa": sigma / 1e6,
        "L_mm_per_end": length * 1e3,
        "submitted_opening_um": explore.delta_a_closed(c, sigma) * 1e6,
        "actual_opening_formula_um": correct_opening(c, sigma) * 1e6,
        "actual_opening_Green_um": reference.opening_by_green(a, c, sigma, a) * 1e6,
        "tip_K_Pa_sqrt_m": reference.tip_by_weight(c, sigma, a),
        "tip_G_J_per_m2": reference.tip_by_weight(c, sigma, a) ** 2 / ep,
        "load_branch_slope_at_submitted_root": slope(u_wrong),
        "minimum_load_branch_slope_up_to_submitted_root": min(map(slope, u_sample)),
    },
    "opening_at_reference_event": {
        "actual_um": correct_opening(
            reference.A + target_length / 1e3, target_sigma * 1e6
        ) * 1e6,
        "submitted_formula_um": explore.delta_a_closed(
            reference.A + target_length / 1e3, target_sigma * 1e6
        ) * 1e6,
    },
    "valid_K_quadrature_check": {
        "closed": K_exact, "quadrature": K_quad,
        "relative_difference": abs(K_quad / K_exact - 1),
    },
    "finite_opening_crosscheck_example": {
        "c_mm": example_c * 1e3,
        "Sigma_MPa": example_sigma / 1e6,
        "submitted_delta_a_um": explore.delta_a_closed(example_c, example_sigma) * 1e6,
        "actual_delta_a_um": correct_opening(example_c, example_sigma) * 1e6,
        "Green_delta_a_um": reference.opening_by_green(a, example_c, example_sigma, a) * 1e6,
        "submitted_full_field_center_gap_um": submitted_center_gap * 1e6,
        "Green_center_gap_um": reference.opening_by_green(0, example_c, example_sigma, a) * 1e6,
    },
    "branch_claim_audit": {
        "stress_maximum_L_mm": a * math.expm1(u_peak) * 1e3,
        "stress_maximum_MPa": stress(a * math.exp(u_peak)) / 1e6,
        "reference_event_u": u_reference,
        "submitted_event_u": u_wrong,
        "claimed_saturation_um": 4 * t * a / ep * 1e6,
        "submitted_formula_zero_Gamma_asymptote_um": 4 * t * a / ep * (1 - 2 / math.pi) * 1e6,
        "submitted_formula_nonzero_Gamma_asymptote": "diverges as 4*K0*sqrt(c)/(sqrt(pi)*E_prime)",
    },
}

print(json.dumps(result, indent=2))
