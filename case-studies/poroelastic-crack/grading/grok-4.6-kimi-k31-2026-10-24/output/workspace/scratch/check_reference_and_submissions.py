#!/usr/bin/env python3
"""Independent numerical checks for the poroelastic-crack grading task.

Does not import or execute the submissions in place. Reconstructs the
corrected-key numbers, the historical paper numbers, and the A/B reported
values from the formulas stated in the submissions.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

NU = 0.2
NU_U = 0.4
ALPHA = 0.89

# Detournay–Cheng plane-strain coupling coefficient (REFERENCE.md Eq. A3).
ETA = ALPHA * (1 - 2 * NU) / (2 * (1 - NU))


def close(actual: float, expected: float, rel: float = 1e-12, abs_: float = 1e-12) -> bool:
    return math.isclose(actual, expected, rel_tol=rel, abs_tol=abs_)


def within_grade_tol(actual: float, expected: float, *, is_k2_0: bool = False) -> bool:
    if is_k2_0:
        return abs(actual - expected) <= 0.01
    if expected == 0.0:
        return abs(actual - expected) <= 0.01
    return abs(actual - expected) <= 0.01 * abs(expected)


# --- Corrected key (REFERENCE.md, P374) ---
corrected = {
    "K1_0_tip_first": 0.807593,
    "K1_0_undrained_first": 1.0,
    "K1_inf": 1.0,
    "V1_0": math.pi * (1 - NU_U),
    "V1_inf": math.pi * (1 - NU),
    # faces-only paper coefficient plus tip contribution 1.086165
    "V1_short_faces": 16 * (1 - NU_U) * (NU_U - NU) / ((1 - NU) * math.sqrt(math.pi)),
    "V1_short_tips": 1.086165,
    "V1_short_coeff": 2.440220,
    "p_faces_0": 2 * (NU_U - NU) / (ALPHA * (1 - 2 * NU)),
    "K2_0": 0.0,
    "K2_inf": -ETA,
    "V2_inf": -ALPHA * math.pi * (1 - 2 * NU) / 2.0,
}

# Historical paper values kept on record, not used as the numerical key.
paper = {
    "K1_0": (1 - NU_U) / (1 - NU),  # 0.75, heuristic
    "V1_short_coeff": 16 * (1 - NU_U) * (NU_U - NU) / ((1 - NU) * math.sqrt(math.pi)),
    "K2_inf": -ETA,
    "V2_inf": -ALPHA * math.pi * (1 - 2 * NU) / 2.0,
}

# Formulas as stated by the submissions (normalized K / sqrt(pi L), G V / L^2).
A_from_formulas = {
    "K1_0": 1.0,
    "K1_inf": 1.0,
    "V1_0": math.pi * (1 - NU_U),
    "V1_inf": math.pi * (1 - NU),
    "V1_short_coeff": 16 * (NU_U - NU) * (1 - NU_U) / ((1 - NU) * math.sqrt(math.pi)),
    "p_faces_0": 2 * (NU_U - NU) / (ALPHA * (1 - 2 * NU)),
    "K2_0": 0.0,
    "K2_inf": -ALPHA,
    "V2_inf": -ALPHA * math.pi * (1 - NU),
}

B_from_formulas = {
    "K1_0": 1.0,
    "K1_inf": 1.0,
    "V1_0": math.pi * (1 - NU_U),
    "V1_inf": math.pi * (1 - NU),
    "V1_short_coeff": 8 * (NU_U - NU) / ((1 - NU) * math.sqrt(math.pi)),
    "p_faces_0": 2 * (NU_U - NU) / (ALPHA * (1 - 2 * NU)),
    "K2_0": 0.0,
    "K2_inf": -ALPHA,
    "V2_inf": -ALPHA * math.pi * (1 - NU),
}

workspace = Path(__file__).resolve().parents[1]
A_json = json.loads((workspace / "A" / "results.json").read_text())
B_json = json.loads((workspace / "B" / "results.json").read_text())

# A's inner/outer split of the faces-only coefficient (A/verify.py lines 137-139).
D = NU_U - NU
A_inner = 8 * D / ((1 - NU) * math.sqrt(math.pi))
A_outer = 8 * D * (1 - 2 * NU_U) / ((1 - NU) * math.sqrt(math.pi))

# Equivalent loading-2 identities.
V2_inf_via_eta = -ETA * math.pi * (1 - NU)
V2_inf_via_alpha_half = -ALPHA * math.pi * (1 - 2 * NU) / 2.0

report: dict = {
    "material": {"nu": NU, "nu_u": NU_U, "alpha": ALPHA, "eta": ETA},
    "corrected_key": corrected,
    "paper_historical": {
        "K1_0": paper["K1_0"],
        "V1_short_coeff": paper["V1_short_coeff"],
        "K2_inf": paper["K2_inf"],
        "V2_inf": paper["V2_inf"],
    },
    "identities": {
        "eta_equals_alpha_times_plane_strain_factor": close(
            ETA, ALPHA * (1 - 2 * NU) / (2 * (1 - NU))
        ),
        "V2_inf_eta_form_matches_alpha_form": close(V2_inf_via_eta, V2_inf_via_alpha_half),
        "V1_inf_over_V1_0": corrected["V1_inf"] / corrected["V1_0"],
        "V1_inf_over_V1_0_expected": (1 - NU) / (1 - NU_U),
        "paper_K1_0_is_075": close(paper["K1_0"], 0.75),
        "faces_plus_tips": corrected["V1_short_faces"] + corrected["V1_short_tips"],
        "A_inner": A_inner,
        "A_outer": A_outer,
        "A_inner_plus_outer": A_inner + A_outer,
        "A_inner_equals_B_short": close(A_inner, B_from_formulas["V1_short_coeff"]),
        "A_faces_equals_paper": close(A_from_formulas["V1_short_coeff"], paper["V1_short_coeff"]),
        "combined_K_A_B": A_from_formulas["K1_inf"] + A_from_formulas["K2_inf"],
        "combined_K_correct": 1.0 + corrected["K2_inf"],
        "combined_V_A_B": A_from_formulas["V1_inf"] + A_from_formulas["V2_inf"],
        "combined_V_correct": corrected["V1_inf"] + corrected["V2_inf"],
        "naive_effective_combined_K": 1.0 - ALPHA,
        "correct_combined_K": 1.0 - ETA,
    },
    "json_matches_stated_formulas": {
        "A": {},
        "B": {},
    },
    "numerical_agreement_vs_corrected": {
        "A": {},
        "B": {},
    },
    "numerical_agreement_vs_paper_historical": {
        "A": {},
        "B": {},
    },
}

numeric_keys = [
    "K1_0",
    "K1_inf",
    "V1_0",
    "V1_inf",
    "V1_short_coeff",
    "p_faces_0",
    "K2_0",
    "K2_inf",
    "V2_inf",
]

for label, js, formulas in (("A", A_json, A_from_formulas), ("B", B_json, B_from_formulas)):
    for key in numeric_keys:
        report["json_matches_stated_formulas"][label][key] = close(js[key], formulas[key], rel=1e-12, abs_=1e-12)

    # Corrected-key comparison. K1_0 is compared both ways.
    report["numerical_agreement_vs_corrected"][label] = {
        "K1_0_vs_tip_first": {
            "reported": js["K1_0"],
            "corrected_reference": corrected["K1_0_tip_first"],
            "agrees": within_grade_tol(js["K1_0"], corrected["K1_0_tip_first"]),
        },
        "K1_0_vs_undrained_first": {
            "reported": js["K1_0"],
            "corrected_reference": corrected["K1_0_undrained_first"],
            "agrees": within_grade_tol(js["K1_0"], corrected["K1_0_undrained_first"]),
        },
        "K1_inf": {
            "reported": js["K1_inf"],
            "corrected_reference": corrected["K1_inf"],
            "agrees": within_grade_tol(js["K1_inf"], corrected["K1_inf"]),
        },
        "V1_0": {
            "reported": js["V1_0"],
            "corrected_reference": corrected["V1_0"],
            "agrees": within_grade_tol(js["V1_0"], corrected["V1_0"]),
        },
        "V1_inf": {
            "reported": js["V1_inf"],
            "corrected_reference": corrected["V1_inf"],
            "agrees": within_grade_tol(js["V1_inf"], corrected["V1_inf"]),
        },
        "V1_short_coeff": {
            "reported": js["V1_short_coeff"],
            "corrected_reference": corrected["V1_short_coeff"],
            "agrees": within_grade_tol(js["V1_short_coeff"], corrected["V1_short_coeff"]),
        },
        "p_faces_0": {
            "reported": js["p_faces_0"],
            "corrected_reference": corrected["p_faces_0"],
            "agrees": within_grade_tol(js["p_faces_0"], corrected["p_faces_0"]),
        },
        "K2_0": {
            "reported": js["K2_0"],
            "corrected_reference": corrected["K2_0"],
            "agrees": within_grade_tol(js["K2_0"], corrected["K2_0"], is_k2_0=True),
        },
        "K2_inf": {
            "reported": js["K2_inf"],
            "corrected_reference": corrected["K2_inf"],
            "agrees": within_grade_tol(js["K2_inf"], corrected["K2_inf"]),
        },
        "V2_inf": {
            "reported": js["V2_inf"],
            "corrected_reference": corrected["V2_inf"],
            "agrees": within_grade_tol(js["V2_inf"], corrected["V2_inf"]),
        },
    }
    report["numerical_agreement_vs_paper_historical"][label] = {
        "K1_0_paper_heuristic": {
            "reported": js["K1_0"],
            "paper": paper["K1_0"],
            "agrees": within_grade_tol(js["K1_0"], paper["K1_0"]),
        },
        "V1_short_coeff_paper_faces_only": {
            "reported": js["V1_short_coeff"],
            "paper": paper["V1_short_coeff"],
            "agrees": within_grade_tol(js["V1_short_coeff"], paper["V1_short_coeff"]),
        },
    }


# Sanity: reconstructed corrected numbers vs REFERENCE.md printed digits.
sanity = {
    "V1_0_printed": 1.88496,
    "V1_inf_printed": 2.51327,
    "p_faces_0_printed": 0.74906,
    "K2_inf_printed": -0.33375,
    "V2_inf_printed": -0.83881,
    "paper_short_printed": 1.35406,
    "eta_printed": 0.33375,
}
sanity_ok = {
    "V1_0": abs(corrected["V1_0"] - sanity["V1_0_printed"]) < 5e-6,
    "V1_inf": abs(corrected["V1_inf"] - sanity["V1_inf_printed"]) < 5e-6,
    "p_faces_0": abs(corrected["p_faces_0"] - sanity["p_faces_0_printed"]) < 5e-6,
    "K2_inf": abs(corrected["K2_inf"] - sanity["K2_inf_printed"]) < 5e-6,
    "V2_inf": abs(corrected["V2_inf"] - sanity["V2_inf_printed"]) < 5e-6,
    "paper_short": abs(paper["V1_short_coeff"] - sanity["paper_short_printed"]) < 5e-6,
    "eta": abs(ETA - sanity["eta_printed"]) < 5e-6,
    "faces_plus_tips_vs_key": abs(
        corrected["V1_short_faces"] + corrected["V1_short_tips"] - corrected["V1_short_coeff"]
    )
    < 5e-7,
}
report["reference_digit_sanity"] = {
    "computed": {
        "V1_0": corrected["V1_0"],
        "V1_inf": corrected["V1_inf"],
        "p_faces_0": corrected["p_faces_0"],
        "K2_inf": corrected["K2_inf"],
        "V2_inf": corrected["V2_inf"],
        "paper_short": paper["V1_short_coeff"],
        "eta": ETA,
        "faces_plus_tips": corrected["V1_short_faces"] + corrected["V1_short_tips"],
    },
    "ok": sanity_ok,
}

out = Path(__file__).with_name("check_output.json")
out.write_text(json.dumps(report, indent=2) + "\n")

print("eta =", ETA)
print("corrected V1_0, V1_inf, p_faces_0 =", corrected["V1_0"], corrected["V1_inf"], corrected["p_faces_0"])
print("corrected K2_inf, V2_inf =", corrected["K2_inf"], corrected["V2_inf"])
print("paper faces-only A =", paper["V1_short_coeff"])
print("A inner, outer, sum =", A_inner, A_outer, A_inner + A_outer)
print("B short =", B_from_formulas["V1_short_coeff"])
print("faces+tips =", corrected["V1_short_faces"] + corrected["V1_short_tips"])
print("sanity_ok =", sanity_ok)
print("A json vs formulas =", report["json_matches_stated_formulas"]["A"])
print("B json vs formulas =", report["json_matches_stated_formulas"]["B"])
print("A vs corrected agrees:")
for k, v in report["numerical_agreement_vs_corrected"]["A"].items():
    print(" ", k, v["agrees"], "reported", v["reported"], "ref", v["corrected_reference"])
print("B vs corrected agrees:")
for k, v in report["numerical_agreement_vs_corrected"]["B"].items():
    print(" ", k, v["agrees"], "reported", v["reported"], "ref", v["corrected_reference"])
print("wrote", out)
