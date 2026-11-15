#!/usr/bin/env python3
"""Compare saved FEM calculations, check them, and write results.json."""

import json
from pathlib import Path

import numpy as np

from solve import estimate, exact_values


def angles(k):
    return np.degrees(np.angle(np.r_[k, np.sum(k)]))


def angle_difference(a, b):
    return (np.asarray(a) - np.asarray(b) + 180) % 360 - 180


def traction_fit(data, material):
    # Average the two radial Gauss locations after separate extrapolations.
    estimates = []
    for xi in (-1 / np.sqrt(3), 1 / np.sqrt(3)):
        sel = ((data[:, 0].real >= 1e-4) & (data[:, 0].real <= 1e-2)
               & np.isclose(data[:, 1].real, xi))
        r = data[sel, 0].real
        design = np.column_stack([np.ones(len(r)), r, r**2])
        columns = slice(2, 4) if material == "upper" else slice(4, 6)
        estimates.append(np.linalg.lstsq(design, data[sel, columns], rcond=None)[0][0])
    return np.mean(estimates, axis=0)


def main():
    exact = exact_values()
    runs = {}
    for name in ("coarse", "medium", "fine", "short"):
        prefix = Path("fem_" + name)
        metadata = json.loads(prefix.with_suffix(".json").read_text())
        with np.load(prefix.with_suffix(".npz")) as data:
            k = estimate(data["cod"])
            j = data["j"]
            selected_j = j[(j[:, 0] > 0.001) & (j[:, 0] < 0.5), 1:]
            runs[name] = {
                "nodes": metadata["nodes"],
                "phases_deg": angles(k).tolist(),
                "Kabs_cod": np.abs(k).tolist(),
                "Kabs_relative_errors": (np.abs(k) / exact["Kabs"] - 1).tolist(),
                "phase_difference_error_deg": float(
                    angles(k)[0] - angles(k)[1] - exact["phase_difference_abs"]
                ),
                "fit_window_change_deg": angle_difference(
                    angles(estimate(data["cod"], 1e-5, 1e-3)), angles(k)
                ).tolist(),
                "J_relative_error_min": (np.min(selected_j, axis=0) / exact["G"][:2] - 1).tolist(),
                "J_relative_error_max": (np.max(selected_j, axis=0) / exact["G"][:2] - 1).tolist(),
            }
            for side in ("upper", "lower"):
                kt = traction_fit(data["tractions"], side)
                runs[name][f"traction_{side}_phases_deg"] = angles(kt).tolist()
                runs[name][f"traction_{side}_Kabs_relative_errors"] = (
                    np.abs(kt) / exact["Kabs"] - 1
                ).tolist()
        assert np.max(np.abs(metadata["force_balance"])) < 1e-12
        assert np.max(np.abs(metadata["moment_balance"])) < 1e-12
        assert metadata["residual_max"] < 1e-5

    fine = runs["fine"]
    assert np.max(np.abs(fine["Kabs_relative_errors"])) < 2e-5
    assert np.max(np.abs(fine["fit_window_change_deg"])) < 0.01
    assert abs(fine["phase_difference_error_deg"]) < 0.005
    for side in ("upper", "lower"):
        assert np.max(np.abs(angle_difference(fine[f"traction_{side}_phases_deg"], fine["phases_deg"]))) < 0.2
        assert np.max(np.abs(fine[f"traction_{side}_Kabs_relative_errors"])) < 0.01
    assert np.max(np.abs(fine["J_relative_error_min"])) < 0.01
    assert np.max(np.abs(fine["J_relative_error_max"])) < 0.01
    mesh_change = angle_difference(runs["fine"]["phases_deg"], runs["medium"]["phases_deg"])
    length_change = angle_difference(runs["short"]["phases_deg"], runs["medium"]["phases_deg"])
    assert np.max(np.abs(mesh_change)) < 0.01
    assert np.max(np.abs(length_change)) < 0.01

    # Use the exact magnitudes and fine-mesh phases to form the combined load.
    psi = fine["phases_deg"][:2]
    k_exact_mag = np.array(exact["Kabs"]) * np.exp(1j * np.radians(psi))
    psi_both = float(np.degrees(np.angle(np.sum(k_exact_mag))))
    results = {
        "eps": exact["eps"],
        "caseA": {
            "G_E2_h_over_P2": exact["G"][0],
            "Kabs_sqrt_h_over_P": exact["Kabs"][0],
            "psi_deg": round(psi[0], 4),
        },
        "caseB": {
            "G_E2_h3_over_M2": exact["G"][1],
            "Kabs_h32_over_M": exact["Kabs"][1],
            "psi_deg": round(psi[1], 4),
        },
        "psi_both_deg": round(psi_both, 4),
    }
    # Verify the exact requested JSON schema, including JSON-number values.
    assert set(results) == {"eps", "caseA", "caseB", "psi_both_deg"}
    assert set(results["caseA"]) == {"G_E2_h_over_P2", "Kabs_sqrt_h_over_P", "psi_deg"}
    assert set(results["caseB"]) == {"G_E2_h3_over_M2", "Kabs_h32_over_M", "psi_deg"}
    reported_numbers = [results["eps"], *results["caseA"].values(),
                        *results["caseB"].values(), results["psi_both_deg"]]
    assert all(type(value) in (float, int) and np.isfinite(value) for value in reported_numbers)
    report = {
        "exact": exact,
        "runs": runs,
        "fine_minus_medium_phases_deg": mesh_change.tolist(),
        "length12_minus_length24_phases_deg": length_change.tolist(),
        "combined_phase_with_exact_magnitudes_deg": psi_both,
        "all_checks_passed": True,
    }
    Path("checks.json").write_text(json.dumps(report, indent=2) + "\n")
    Path("results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print("Results:")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
